"""How much does an undisturbed fish vary on its own?

This is the number that decides whether behavioural monitoring can work, and
it is almost never reported.  A detector fires when a behavioural metric
leaves its normal range; the size of that normal range, measured on animals
that nothing happened to, is the denominator of every later claim.

Two questions here.

1. **Whose variance is it?**  If a metric differs mostly between individuals
   and hardly at all between recordings of the same individual, then it is a
   trait, and a threshold built on a population will spend its budget on
   telling fish apart rather than on detecting events.  The one-way
   random-effects intraclass correlation gives the split directly:

       ICC = var_between_fish / (var_between_fish + var_within_fish)

2. **Is the distribution anything like Gaussian?**  Behavioural quantities
   are famously not, and a threshold set at "mean plus three standard
   deviations" assumes they are.  Measured here as the excess kurtosis of each
   metric and, more usefully, as the realised exceedance rate of a nominal
   three-sigma rule.  Every factor by which that exceeds 0.27% is a factor of
   extra crew callouts.

Data: 504 usable recordings of spontaneous swimming from 39 larval zebrafish,
after the tracking-quality gate.
"""

import json
import sys

sys.path.insert(0, 'src')

import numpy as np

from behaviour import GROUPS, METRICS, load, state_matrix, valid

SIGMAS = [2.0, 3.0, 4.0]


def icc_oneway(values, groups):
    """One-way random-effects ICC(1) with unequal group sizes."""
    g = np.unique(groups)
    k = len(g)
    n_i = np.array([np.sum(groups == x) for x in g], float)
    N = n_i.sum()
    means = np.array([values[groups == x].mean() for x in g])
    grand = values.mean()
    ss_between = float(np.sum(n_i * (means - grand) ** 2))
    ss_within = float(np.sum([np.sum((values[groups == x]
                                      - values[groups == x].mean()) ** 2)
                              for x in g]))
    ms_b = ss_between / (k - 1)
    ms_w = ss_within / (N - k)
    n0 = (N - np.sum(n_i ** 2) / N) / (k - 1)      # effective group size
    var_b = max((ms_b - ms_w) / n0, 0.0)
    icc = var_b / (var_b + ms_w) if (var_b + ms_w) > 0 else 0.0
    return {"icc": float(icc), "var_between": float(var_b),
            "var_within": float(ms_w), "n_groups": int(k), "n": int(N)}


def tail_report(v, sigmas=SIGMAS):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    m, s = v.mean(), v.std(ddof=1)
    from math import erfc, sqrt
    out = {"excess_kurtosis": float(((v - m) ** 4).mean() / s ** 4 - 3.0)}
    for k in sigmas:
        nominal = erfc(k / sqrt(2.0))            # two-sided
        realised = float(np.mean(np.abs(v - m) > k * s))
        out[f"{k:g}sigma"] = {
            "nominal": nominal, "realised": realised,
            "ratio": realised / nominal if nominal > 0 else np.nan}
    return out


def main():
    trials = [t for t in load() if valid(t)]
    X, fish, names = state_matrix(trials)
    dur = np.array([t["t"].max() - t["t"].min() for t in trials])

    out = {"config": {
        "n_recordings": int(len(trials)),
        "n_fish": int(len(np.unique(fish))),
        "n_bouts": int(sum(len(t["t"]) for t in trials)),
        "median_duration_s": float(np.median(dur)),
        "total_hours": float(dur.sum() / 3600.0),
        "sigmas": SIGMAS,
    }, "icc": {}, "tails": {}}

    for j, m in enumerate(names):
        v = X[:, j]
        ok = np.isfinite(v)
        out["icc"][m] = icc_oneway(v[ok], fish[ok])
        out["tails"][m] = tail_report(v[ok])
        out["icc"][m]["group"] = next(g for g, ms in GROUPS.items() if m in ms)

    # the same two questions for the raw bout-level quantities, which is what a
    # frame-by-frame detector would threshold
    raw = {"turn_angle_deg": np.concatenate([t["turn"] for t in trials]),
           "displacement": np.abs(np.concatenate([t["r"] for t in trials])),
           "inter_bout_interval_s": np.concatenate(
               [np.diff(t["t"]) for t in trials])}
    out["tails_bout_level"] = {k: tail_report(v) for k, v in raw.items()}

    with open("results/exp1_variability.json", "w") as fh:
        json.dump(out, fh, indent=1)

    c = out["config"]
    print(f"{c['n_recordings']} recordings, {c['n_fish']} fish, "
          f"{c['n_bouts']} bouts, {c['total_hours']:.1f} h of behaviour\n")
    print(f"{'metric':24s} {'group':10s} {'ICC':>6s} "
          f"{'exc.kurt':>9s} {'3s real':>8s} {'x nominal':>10s}")
    for m in names:
        i, t = out["icc"][m], out["tails"][m]
        print(f"{m:24s} {i['group']:10s} {i['icc']:6.2f} "
              f"{t['excess_kurtosis']:9.1f} {t['3sigma']['realised']*100:7.2f}% "
              f"{t['3sigma']['ratio']:10.1f}")
    print("\nbout-level quantities (what a per-event threshold would see)")
    for k, t in out["tails_bout_level"].items():
        print(f"  {k:24s} excess kurtosis {t['excess_kurtosis']:8.1f}   "
              f"3-sigma exceedance {t['3sigma']['realised']*100:5.2f}% "
              f"= {t['3sigma']['ratio']:.0f}x nominal")


if __name__ == "__main__":
    main()
