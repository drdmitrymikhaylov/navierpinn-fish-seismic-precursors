"""Checks on the loader, the metrics and the estimators."""

import sys

sys.path.insert(0, 'src')

import numpy as np

from behaviour import (METRICS, circular_dispersion, load, state,
                       state_matrix, valid)
from exp1_variability import icc_oneway, tail_report
from exp2_detector import robust_scale, stat_scan, stat_single


def test_dataset_shape():
    trials = load()
    assert len(trials) > 400
    good = [t for t in trials if valid(t)]
    assert 400 < len(good) < len(trials)
    fish = {t["fish"] for t in good}
    assert len(fish) > 30
    for t in good:
        assert len(t["t"]) >= 8
        assert np.all(np.diff(t["t"]) > 0), "bout times must increase"
        assert len(t["turn"]) >= 7


def test_quality_gate_rejects_out_of_arena():
    """A track outside the arena is a tracking failure, not a fish."""
    trials = [t for t in load() if valid(t)]
    bad = dict(trials[0])
    bad["x"] = bad["x"].copy()
    bad["x"][0] = 5000.0
    assert not valid(bad)


def test_circular_dispersion_bounds():
    assert circular_dispersion(np.zeros(50)) < 1e-12
    rng = np.random.default_rng(0)
    uniform = circular_dispersion(rng.uniform(-180, 180, 200000))
    assert 0.99 < 1 - abs(1 - uniform) <= 1.0
    assert uniform > 0.99


def test_metrics_are_finite_and_sane():
    trials = [t for t in load() if valid(t)]
    X, fish, names = state_matrix(trials)
    assert names == METRICS
    assert np.isfinite(X).all()
    b = X[:, names.index("bout_rate")]
    assert 0.1 < np.median(b) < 5.0, np.median(b)
    p = X[:, names.index("turn_persistence")]
    assert 0.0 <= p.min() and p.max() <= 1.0
    f = X[:, names.index("large_turn_fraction")]
    assert 0.0 <= f.min() and f.max() <= 1.0


def test_icc_recovers_a_known_split():
    """Synthetic data with a known between-group share must return it."""
    rng = np.random.default_rng(0)
    k, n = 40, 15
    for true_icc in (0.2, 0.6):
        sb = np.sqrt(true_icc)
        sw = np.sqrt(1 - true_icc)
        groups = np.repeat(np.arange(k), n)
        vals = np.repeat(rng.normal(0, sb, k), n) + rng.normal(0, sw, k * n)
        est = icc_oneway(vals, groups)["icc"]
        assert abs(est - true_icc) < 0.1, (true_icc, est)


def test_tail_report_is_calibrated_on_gaussian_data():
    """The exceedance ratio must be about 1 when the data really are normal."""
    v = np.random.default_rng(1).normal(0, 1, 400000)
    r = tail_report(v)
    assert abs(r["excess_kurtosis"]) < 0.1
    assert 0.7 < r["3sigma"]["ratio"] < 1.4, r["3sigma"]


def test_behaviour_really_is_heavier_tailed_than_that():
    """The central claim of the repository, pinned."""
    trials = [t for t in load() if valid(t)]
    turn = np.concatenate([t["turn"] for t in trials])
    ibi = np.concatenate([np.diff(t["t"]) for t in trials])
    assert tail_report(turn)["excess_kurtosis"] > 20
    assert tail_report(ibi)["excess_kurtosis"] > 50
    X, fish, names = state_matrix(trials)
    ratios = [tail_report(X[:, j])["3sigma"]["ratio"] for j in range(X.shape[1])]
    assert max(ratios) > 5.0, max(ratios)


def test_single_statistic_is_never_above_the_scan():
    trials = [t for t in load() if valid(t)]
    X, fish, names = state_matrix(trials)
    med, mad = robust_scale(X)
    j = names.index("bout_rate")
    assert np.all(stat_single(X, med, mad, j) <= stat_scan(X, med, mad) + 1e-12)


def test_drift_diffusion_matches_binned_moments():
    try:
        import torch  # noqa: F401
    except ImportError:
        print("  (skipped: no torch)"); return
    from exp3_fp_pinn import (binned_moments, fit_drift_diffusion, stationary,
                              transitions)
    trials = [t for t in load() if valid(t)]
    r, dr, dt = transitions(trials)
    m, _ = fit_drift_diffusion(r, dr, dt, steps=2000)
    b = binned_moments(r, dr, dt)
    import torch
    with torch.no_grad():
        R = torch.tensor([[x["r"]] for x in b])
        a_fit = m.a(R).squeeze().numpy(); D_fit = m.D(R).squeeze().numpy()
    # the drift is not monotonic in r, so Pearson on the bin values, not rank
    assert np.corrcoef([x["a"] for x in b], a_fit)[0, 1] > 0.85
    assert np.corrcoef([x["D"] for x in b], D_fit)[0, 1] > 0.85
    rg, p, _, _ = stationary(m)
    assert abs(np.trapezoid(p, rg) - 1) < 1e-6


def test_logratio_pinn_matches_finite_differences():
    try:
        import torch  # noqa: F401
    except ImportError:
        print("  (skipped: no torch)"); return
    import torch
    from exp3_fp_pinn import (fd_solution, fit_drift_diffusion, pinn_solution,
                              stationary, transitions)
    trials = [t for t in load() if valid(t)]
    r, dr, dt = transitions(trials)
    m, _ = fit_drift_diffusion(r, dr, dt, steps=800)
    rg, p_stat, _, _ = stationary(m)
    p0 = lambda x: np.exp(-0.5 * ((x - 0.2) / 0.08) ** 2) + 1e-3
    times = np.linspace(0, 300.0, 7)
    rc, P_fd = fd_solution(m, p0, times)
    p_of, _ = pinn_solution(m, p0, (rg, np.log(p_stat + 1e-300)), steps=2500)
    with torch.no_grad():
        RR, TT = np.meshgrid(rc, times, indexing="ij")
        P = p_of(torch.tensor(RR.reshape(-1, 1)), torch.tensor(TT.reshape(-1, 1))).numpy().reshape(RR.shape)
    assert np.max(np.abs(P - P_fd)) / np.max(P_fd) < 0.05


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
