"""Space use as a physical process: drift, diffusion and a Fokker-Planck PINN.

Thigmotaxis -- the tendency to stay near the wall -- is the standard
behavioural readout for anxiety-like states in fish, and it is usually
reported as one number, the median distance from the centre.  That number
is the stationary state of a physical process: a fish is a random walker
in a bounded arena with a position-dependent drift toward the wall and a
position-dependent diffusion.  Written that way, the readout has dynamics,
and the dynamics say how long a change takes to become visible.

Three steps, each checked against something that does not use a network:

  1. learn the drift a(r) and diffusion D(r) of the radial coordinate from
     the bout-to-bout transitions (neural Kramers-Moyal: two small networks
     fitted by the Gaussian transition likelihood), and compare them with
     the binned conditional moments

  2. the stationary Fokker-Planck solution predicted by a(r), D(r) is
     compared with the observed distribution of radius: does the process
     reproduce thigmotaxis it was not fitted to?  (It was fitted to
     transitions, not to where the fish are.)

  3. a PINN solves the time-dependent Fokker-Planck equation
        dp/dt = -d/dr [a(r) p] + d^2/dr^2 [D(r) p]
     with reflecting walls, from a fish released near the centre, and is
     checked against a finite-difference solution.  What the solution buys
     is a relaxation time: how many seconds until the radial distribution
     has forgotten where the fish started.  That is the shortest window in
     which a change in space use can in principle be seen, and it is
     compared with the seven minutes exp2 found empirically.

Units: radius in arena half-widths (0 = centre, 1 = wall), time in seconds.
"""

import json
import sys

sys.path.insert(0, 'src')

import numpy as np
import torch
from scipy.integrate import solve_ivp

from behaviour import load, valid

torch.set_default_dtype(torch.float64)
torch.manual_seed(0)
R_MAX = 1.0
T_MAX = 300.0


def transitions(trials):
    r, dr, dt = [], [], []
    for t in trials:
        rr = np.hypot(t["x"] - 500.0, t["y"] - 500.0) / 500.0
        r.append(rr[:-1]); dr.append(np.diff(rr)); dt.append(np.diff(t["t"]))
    r, dr, dt = (np.concatenate(v) for v in (r, dr, dt))
    keep = (dt > 0.05) & (dt < 20) & (r < R_MAX)
    return r[keep], dr[keep], dt[keep]


def mlp(width=32, depth=2, n_out=1):
    layers, d = [], 1
    for _ in range(depth):
        layers += [torch.nn.Linear(d, width), torch.nn.Tanh()]
        d = width
    layers.append(torch.nn.Linear(d, n_out))
    return torch.nn.Sequential(*layers)


class DriftDiffusion(torch.nn.Module):
    """a(r) free; D(r) = softplus(.) > 0."""

    def __init__(self):
        super().__init__()
        self.a_net, self.d_net = mlp(), mlp()

    def a(self, r):
        return 0.05 * self.a_net(r)

    def D(self, r):
        return 2e-3 * torch.nn.functional.softplus(self.d_net(r) + 1.0)


def fit_drift_diffusion(r, dr, dt, steps=3000):
    m = DriftDiffusion()
    R = torch.tensor(r).reshape(-1, 1)
    DR = torch.tensor(dr).reshape(-1, 1)
    DT = torch.tensor(dt).reshape(-1, 1)
    opt = torch.optim.Adam(m.parameters(), lr=3e-3)
    for _ in range(steps):
        opt.zero_grad()
        mu = m.a(R) * DT
        var = 2.0 * m.D(R) * DT
        nll = torch.mean(0.5 * torch.log(var) + (DR - mu) ** 2 / (2 * var))
        nll.backward()
        opt.step()
    return m, float(nll)


def binned_moments(r, dr, dt, nb=10):
    edges = np.linspace(0, R_MAX, nb + 1)
    out = []
    for i in range(nb):
        k = (r >= edges[i]) & (r < edges[i + 1])
        if k.sum() < 30:
            continue
        out.append({"r": float(0.5 * (edges[i] + edges[i + 1])), "n": int(k.sum()),
                    "a": float(np.mean(dr[k] / dt[k])),
                    "D": float(np.mean(dr[k] ** 2 / dt[k]) / 2)})
    return out


def stationary(m, n=400):
    """Zero-flux stationary solution: p ~ (1/D) exp(int a/D dr)."""
    r = np.linspace(1e-3, R_MAX, n)
    with torch.no_grad():
        R = torch.tensor(r).reshape(-1, 1)
        a, D = m.a(R).squeeze().numpy(), m.D(R).squeeze().numpy()
    phi = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] / D[1:] + a[:-1] / D[:-1]) * np.diff(r))])
    p = np.exp(phi - phi.max()) / D
    p /= np.trapezoid(p, r)
    return r, p, a, D


def fd_solution(m, p0_fn, times, n=200):
    """Method-of-lines finite differences, conservative form, reflecting walls."""
    r = np.linspace(0, R_MAX, n + 1)
    rc = 0.5 * (r[1:] + r[:-1]); h = r[1] - r[0]
    with torch.no_grad():
        a_f = m.a(torch.tensor(r).reshape(-1, 1)).squeeze().numpy()      # at faces
        D_c = m.D(torch.tensor(rc).reshape(-1, 1)).squeeze().numpy()     # at centres
    p0 = p0_fn(rc); p0 /= np.sum(p0) * h

    def rhs(t, p):
        # flux at interior faces: J = a p - d(Dp)/dr ; zero at both walls
        pf = 0.5 * (p[1:] + p[:-1])
        dDp = (D_c[1:] * p[1:] - D_c[:-1] * p[:-1]) / h
        J = a_f[1:-1] * pf - dDp
        J = np.concatenate([[0.0], J, [0.0]])
        return -(J[1:] - J[:-1]) / h
    sol = solve_ivp(rhs, (0, times[-1]), p0, t_eval=times, method="BDF",
                    rtol=1e-6, atol=1e-9)
    return rc, sol.y


def pinn_solution(m, p0_fn, log_ps_grid, steps=6000):
    """PINN for the log-ratio u = ln(p / p_stat).

    A network for p itself failed twice.  First it found the trivial
    solution p = 0 for t > 0 (a soft initial condition and sparse
    collocation near t = 0 let it).  With a hard initial condition and a
    mass-conservation loss it then failed to build the boundary layer at
    the wall, where the stationary density peaks at 14 within 5% of the
    arena -- the usual spectral bias of a smooth network.

    The fix is the change of variable.  With p = p_s e^u and p_s the
    stationary solution (a p_s = d(D p_s)/dr), the Fokker-Planck equation
    becomes
        u_t = D u_rr + a u_r + D u_r^2,
    which contains no p_s at all, has reflecting walls as u_r = 0, and whose
    solution is smooth and tends to u = 0.  The boundary layer lives in p_s,
    which is known in closed form; the network only has to learn the
    smooth part.
    """
    rg, lps = log_ps_grid
    net = torch.nn.Sequential(torch.nn.Linear(2, 64), torch.nn.Tanh(),
                              torch.nn.Linear(64, 64), torch.nn.Tanh(),
                              torch.nn.Linear(64, 64), torch.nn.Tanh(),
                              torch.nn.Linear(64, 1))

    def u_of(r, t_s):
        return net(torch.cat([r, t_s], 1))

    r_c = torch.rand(4000, 1) * R_MAX; t_c = torch.rand(4000, 1)
    r_c.requires_grad_(True); t_c.requires_grad_(True)
    # initial condition on a grid (soft): u0 = ln p0 - ln p_s
    r_i = torch.tensor(rg).reshape(-1, 1)
    p0 = p0_fn(rg); p0 = p0 / np.trapezoid(p0, rg)
    u_i = torch.tensor(np.log(p0 + 1e-12) - lps).reshape(-1, 1)
    t_i = torch.zeros_like(r_i)
    ps_i = torch.tensor(np.exp(lps)).reshape(-1, 1)
    t_b = torch.rand(300, 1)
    r_b0 = torch.zeros_like(t_b, requires_grad=True); r_b1 = torch.full_like(t_b, R_MAX, requires_grad=True)
    t_m = torch.rand(40, 1)
    for q in m.parameters():
        q.requires_grad_(False)

    def losses():
        u = u_of(r_c, t_c)
        u_r = torch.autograd.grad(u, r_c, torch.ones_like(u), create_graph=True)[0]
        u_rr = torch.autograd.grad(u_r, r_c, torch.ones_like(u_r), create_graph=True)[0]
        u_t = torch.autograd.grad(u, t_c, torch.ones_like(u), create_graph=True)[0]
        res = u_t - T_MAX * (m.D(r_c) * u_rr + m.a(r_c) * u_r + m.D(r_c) * u_r ** 2)
        l_pde = torch.mean(res ** 2)
        l_ic = torch.mean((u_of(r_i, t_i) - u_i) ** 2)
        # reflecting walls: u_r = 0
        ub0 = u_of(r_b0, t_b); ub1 = u_of(r_b1, t_b)
        g0 = torch.autograd.grad(ub0, r_b0, torch.ones_like(ub0), create_graph=True)[0]
        g1 = torch.autograd.grad(ub1, r_b1, torch.ones_like(ub1), create_graph=True)[0]
        l_bc = torch.mean(g0 ** 2) + torch.mean(g1 ** 2)
        # mass conservation of p = p_s e^u
        pm = torch.stack([(ps_i * torch.exp(u_of(r_i, torch.full_like(r_i, float(tm))))).squeeze()
                          for tm in t_m])
        mass = torch.trapezoid(pm, r_i.squeeze(), dim=1)
        l_mass = torch.mean((mass - 1.0) ** 2)
        return l_pde, l_ic, l_bc, l_mass

    def total():
        a, b, c, d = losses()
        return a + 10 * b + 10 * c + 10 * d

    opt = torch.optim.Adam(net.parameters(), lr=2e-3)
    for it in range(steps):
        opt.zero_grad()
        l = total()
        l.backward()
        opt.step()
    lb = torch.optim.LBFGS(net.parameters(), max_iter=1500, line_search_fn="strong_wolfe")

    def closure():
        lb.zero_grad()
        l = total()
        l.backward()
        return l
    lb.step(closure)
    a, b, c, d = losses()

    def p_of(r, t):
        with torch.no_grad():
            lp = torch.tensor(np.interp(r.squeeze().numpy(), rg, lps)).reshape(-1, 1)
            return torch.exp(lp + u_of(r, t / T_MAX))
    return p_of, {"pde": float(a), "ic": float(b), "bc": float(c), "mass": float(d)}


def relaxation_time(rc, P, times, p_stat):
    """First time at which the L1 distance to the stationary density is
    below 0.1 (10% of probability mass misplaced)."""
    h = rc[1] - rc[0]
    d = [0.5 * np.sum(np.abs(P[:, k] / (np.sum(P[:, k]) * h) - p_stat)) * h for k in range(len(times))]
    d = np.array(d)
    below = np.where(d < 0.1)[0]
    return (float(times[below[0]]) if len(below) else None), d.tolist()


def main():
    trials = [t for t in load() if valid(t)]
    r, dr, dt = transitions(trials)
    m, nll = fit_drift_diffusion(r, dr, dt)
    binned = binned_moments(r, dr, dt)
    rg, p_stat, a_g, D_g = stationary(m)

    # observed radial distribution (all bouts, all fish)
    r_all = np.concatenate([np.hypot(t["x"] - 500, t["y"] - 500) / 500 for t in trials])
    hist, edges = np.histogram(r_all, bins=25, range=(0, 1), density=True)
    cent = 0.5 * (edges[1:] + edges[:-1])
    p_model = np.interp(cent, rg, p_stat)
    kl = float(np.sum(hist[hist > 0] * np.log(hist[hist > 0] / p_model[hist > 0])) * (edges[1] - edges[0]))
    med_obs = float(np.median(r_all))
    cdf = np.cumsum(p_stat) * (rg[1] - rg[0]); med_model = float(rg[np.searchsorted(cdf, 0.5)])

    # time-dependent: released near the centre
    # a small floor keeps ln p0 finite for the log-ratio network
    p0 = lambda x: np.exp(-0.5 * ((x - 0.2) / 0.08) ** 2) + 1e-3
    times = np.linspace(0, T_MAX, 61)
    rc, P_fd = fd_solution(m, p0, times)
    p_stat_c = np.interp(rc, rg, p_stat)
    tau_fd, dist_fd = relaxation_time(rc, P_fd, times, p_stat_c)

    p_of, lres = pinn_solution(m, p0, (rg, np.log(p_stat + 1e-300)))
    with torch.no_grad():
        RR, TT = np.meshgrid(rc, times, indexing="ij")
        P_pinn = p_of(torch.tensor(RR.reshape(-1, 1)), torch.tensor(TT.reshape(-1, 1))).numpy().reshape(RR.shape)
    err = float(np.max(np.abs(P_pinn - P_fd)) / np.max(P_fd))
    tau_pinn, dist_pinn = relaxation_time(rc, P_pinn, times, p_stat_c)

    # where do recordings start, and where do they end?  The stationary
    # state is only meaningful if recordings are long enough to reach it.
    r_first = np.array([np.hypot(t["x"][0] - 500, t["y"][0] - 500) / 500 for t in trials])
    r_last = np.array([np.hypot(t["x"][-1] - 500, t["y"][-1] - 500) / 500 for t in trials])
    dur = np.array([t["t"][-1] - t["t"][0] for t in trials])

    # a stressed fish: drift toward the wall 50% stronger.  Starting from the
    # observed first-bout distribution, how far apart are the median radii
    # after one recording (26 s), in units of the between-recording spread?
    hist0, e0 = np.histogram(r_first, bins=20, range=(0, 1), density=True)
    c0 = 0.5 * (e0[1:] + e0[:-1])
    p_init = lambda x: np.interp(x, c0, hist0) + 1e-3
    t26 = np.linspace(0, 26.0, 27)

    class Scaled(torch.nn.Module):
        def __init__(self, f): super().__init__(); self.f = f
        def a(self, rr): return self.f * m.a(rr)
        def D(self, rr): return m.D(rr)

    med26 = {}
    for f in (1.0, 1.5, 2.0):
        rcx, Px = fd_solution(Scaled(f), p_init, t26)
        cdfx = np.cumsum(Px, axis=0) * (rcx[1] - rcx[0])
        med26[str(f)] = [float(rcx[np.searchsorted(cdfx[:, k] / cdfx[-1, k], 0.5)]) for k in range(len(t26))]
    cd = np.array([np.median(np.hypot(t["x"] - 500, t["y"] - 500) / 500) for t in trials])
    robust_sd = float(1.4826 * np.median(np.abs(cd - np.median(cd))))
    effect = {f: float((med26[f][-1] - med26["1.0"][-1]) / robust_sd) for f in ("1.5", "2.0")}
    t_half = None
    out = {"n_transitions": int(len(r)), "nll": nll, "binned": binned,
           "drift_diffusion_curve": {"r": rg.tolist(), "a": a_g.tolist(), "D": D_g.tolist()},
           "stationary": {"r": rg.tolist(), "p": p_stat.tolist(),
                          "observed_centres": cent.tolist(), "observed_density": hist.tolist(),
                          "kl_obs_vs_model": kl, "median_r_observed": med_obs,
                          "median_r_model": med_model},
           "relaxation": {"times": times.tolist(), "tau_fd_s": tau_fd, "tau_pinn_s": tau_pinn,
                          "dist_fd": dist_fd, "dist_pinn": dist_pinn,
                          "pinn_vs_fd_max_rel_err": err, "pinn_losses": lres,
                          "rc": rc.tolist(), "P_fd": P_fd[:, ::10].tolist(),
                          "P_pinn": P_pinn[:, ::10].tolist()},
           "recordings": {"median_r_first_bout": float(np.median(r_first)),
                          "median_r_last_bout": float(np.median(r_last)),
                          "frac_last_bout_at_wall": float(np.mean(r_last > 0.8)),
                          "median_duration_s": float(np.median(dur))},
           "stress_response": {"t_s": t26.tolist(), "median_r_t": med26,
                               "robust_sd_centre_distance": robust_sd,
                               "effect_size_sd_at_26s": effect}}
    with open("results/exp3_fp_pinn.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"{len(r)} transitions; transition NLL {nll:.3f}")
    print("binned  r     a(r)      D(r)       | fitted a, D")
    for b in binned:
        with torch.no_grad():
            R = torch.tensor([[b['r']]])
            print(f"  {b['r']:.2f}  {b['a']:+.4f}  {b['D']:.5f}   | {float(m.a(R)):+.4f}  {float(m.D(R)):.5f}")
    print(f"stationary: KL(observed || model) = {kl:.3f}; median r observed {med_obs:.3f}, "
          f"model {med_model:.3f}")
    print(f"relaxation from centre release: FD {tau_fd} s, PINN {tau_pinn} s; "
          f"PINN vs FD max rel err {err:.3f}; losses {lres}")
    print(f"recordings: first-bout median r {np.median(r_first):.3f}, last-bout "
          f"{np.median(r_last):.3f}, {np.mean(r_last > 0.8)*100:.0f}% end at the wall")
    print(f"after 26 s from the observed start: median r baseline {med26['1.0'][-1]:.3f}, "
          f"drift x1.5 {med26['1.5'][-1]:.3f}, x2 {med26['2.0'][-1]:.3f}; "
          f"effect sizes {effect} (robust SD {robust_sd:.3f})")


if __name__ == "__main__":
    main()
