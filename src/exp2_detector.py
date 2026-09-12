"""A frozen anomaly detector, its real false-alarm rate, and its power.

The detector is deliberately simple, because a detector that will be
pre-registered and published before the data are collected has to be simple
enough to state in one paragraph:

  robust z    each metric is centred on the median and scaled by the median
              absolute deviation of the training recordings
  statistic   two versions, and the difference between them is a result:
              SINGLE, the robust z of one metric fixed in advance, and SCAN,
              the largest absolute robust z across all eleven
  threshold   either the Gaussian rule (three sigma, Bonferroni-corrected for
              the scan) or the empirical (1 - alpha) quantile of the statistic
              on the training data

Fish are split, never recordings: the training fish and the test fish are
disjoint, so nothing the detector learned about one animal helps it on
another.  That is the situation of a deployed system meeting a new tank.

Three things are measured.

  false alarms   how often the detector fires on held-out recordings of
                 undisturbed fish.  The nominal rate is what the rule
                 promises; the realised rate is what it delivers.
  power          a known shift is injected into held-out recordings and the
                 detection rate is measured against the size of that shift.
  observation    the same, after averaging k consecutive recordings, which is
                 what watching for longer actually buys.
"""

import json
import sys

sys.path.insert(0, 'src')

import numpy as np

from behaviour import METRICS, load, state_matrix, valid

ALPHA = 0.01                      # one false alarm in a hundred windows
EFFECTS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]      # in robust SD units
POOL = [1, 2, 4, 8, 16]
N_SPLITS = 200
PRIMARY = "bout_rate"     # the single pre-specified endpoint
RNG = np.random.default_rng(0)


def robust_scale(X):
    med = np.median(X, axis=0)
    mad = np.median(np.abs(X - med), axis=0) * 1.4826
    mad = np.where(mad > 0, mad, np.std(X, axis=0) + 1e-12)
    return med, mad


def zscores(X, med, mad):
    return (X - med) / mad


def stat_scan(X, med, mad):
    """Largest absolute robust z over every metric -- the usual dashboard."""
    return np.max(np.abs(zscores(X, med, mad)), axis=1)


def stat_single(X, med, mad, j):
    """Absolute robust z of one metric chosen before the data were seen."""
    return np.abs(zscores(X, med, mad)[:, j])


def main():
    trials = [t for t in load() if valid(t)]
    X, fish, names = state_matrix(trials)
    dur = np.array([t["t"].max() - t["t"].min() for t in trials])
    keep = np.isfinite(X).all(axis=1)
    X, fish, dur = X[keep], fish[keep], dur[keep]
    fish_ids = np.unique(fish)

    out = {"config": {"alpha": ALPHA, "effects_sd": EFFECTS, "pool": POOL,
                      "n_splits": N_SPLITS, "metrics": names,
                      "n_recordings": int(len(X)), "n_fish": int(len(fish_ids)),
                      "median_duration_s": float(np.median(dur))},
           "false_alarms": {}, "power": {}, "pooling": {}}

    jp = names.index(PRIMARY)
    fa = {"scan|gauss": [], "scan|emp": [], "single|gauss": [], "single|emp": []}
    thr = {"scan": [], "single": []}
    power = {f"{s}|{e}": [] for s in ("scan", "single") for e in EFFECTS}
    pooled = {f"{s}|{k}|{e}": [] for s in ("scan", "single")
              for k in POOL for e in EFFECTS}

    for split in range(N_SPLITS):
        rng = np.random.default_rng(split)
        perm = rng.permutation(fish_ids)
        tr_fish = set(perm[: int(0.6 * len(perm))].tolist())
        tr = np.array([f in tr_fish for f in fish])
        if tr.sum() < 20 or (~tr).sum() < 20:
            continue
        med, mad = robust_scale(X[tr])
        Xte = X[~tr]
        t_emp = {}
        for sname, fn in (("scan", lambda M: stat_scan(M, med, mad)),
                          ("single", lambda M: stat_single(M, med, mad, jp))):
            s_tr, s_te = fn(X[tr]), fn(Xte)
            fa[f"{sname}|gauss"].append(float(np.mean(s_te > 3.0)))
            t = float(np.quantile(s_tr, 1 - ALPHA))
            t_emp[sname] = t
            thr[sname].append(t)
            fa[f"{sname}|emp"].append(float(np.mean(s_te > t)))

        # --- power: shift the primary metric of the held-out recordings.
        # The scan is given the same shift, so the two differ only in how many
        # places they look for it.
        for e in EFFECTS:
            Z = Xte.copy()
            Z[:, jp] = Z[:, jp] + e * mad[jp]
            power[f"scan|{e}"].append(
                float(np.mean(stat_scan(Z, med, mad) > t_emp["scan"])))
            power[f"single|{e}"].append(
                float(np.mean(stat_single(Z, med, mad, jp) > t_emp["single"])))

        # --- pooling: average k held-out recordings before testing.  The
        # threshold has to be recalibrated for each k, because the spread of a
        # mean of k is not the spread of a single record.
        for k in POOL:
            # many random groupings, so the calibration quantile is not itself
            # estimated from a handful of pooled records
            n_draw = 400
            gt = rng.integers(0, tr.sum(), size=(n_draw, k))
            base_tr = X[tr][gt].mean(axis=1)
            gte = rng.integers(0, len(Xte), size=(n_draw, k))
            base = Xte[gte].mean(axis=1)
            med_k, mad_k = robust_scale(base_tr)
            t_k = {"scan": float(np.quantile(
                       stat_scan(base_tr, med_k, mad_k), 1 - ALPHA)),
                   "single": float(np.quantile(
                       stat_single(base_tr, med_k, mad_k, jp), 1 - ALPHA))}
            for e in EFFECTS:
                Z = base.copy()
                Z[:, jp] = Z[:, jp] + e * mad[jp]
                pooled[f"scan|{k}|{e}"].append(
                    float(np.mean(stat_scan(Z, med_k, mad_k) > t_k["scan"])))
                pooled[f"single|{k}|{e}"].append(
                    float(np.mean(stat_single(Z, med_k, mad_k, jp)
                                  > t_k["single"])))

    nominal_gauss = {"single": 0.0027,
                     "scan": 1 - (1 - 0.0027) ** len(names)}
    out["false_alarms"] = {"n_splits": int(len(fa["scan|gauss"])),
                           "primary_metric": PRIMARY}
    for sname in ("single", "scan"):
        out["false_alarms"][sname] = {
            "nominal_gaussian": nominal_gauss[sname],
            "realised_gaussian": float(np.mean(fa[f"{sname}|gauss"])),
            "gaussian_ratio": float(np.mean(fa[f"{sname}|gauss"])
                                    / nominal_gauss[sname]),
            "nominal_empirical": ALPHA,
            "realised_empirical": float(np.mean(fa[f"{sname}|emp"])),
            "empirical_ratio": float(np.mean(fa[f"{sname}|emp"]) / ALPHA),
            "threshold_median": float(np.median(thr[sname])),
        }
    out["power"] = {k: float(np.mean(v)) for k, v in power.items() if v}
    out["pooling"] = {k: float(np.mean(v)) for k, v in pooled.items() if v}

    with open("results/exp2_detector.json", "w") as fh:
        json.dump(out, fh, indent=1)

    f = out["false_alarms"]
    print(f"held-out fish, {f['n_splits']} random splits, "
          f"primary metric = {f['primary_metric']}\n")
    for sname in ("single", "scan"):
        v = f[sname]
        print(f"  [{sname}]")
        print(f"    three-sigma rule   promises "
              f"{v['nominal_gaussian']*100:5.2f}% of windows, delivers "
              f"{v['realised_gaussian']*100:5.2f}%  ({v['gaussian_ratio']:5.1f}x)")
        print(f"    empirical quantile promises "
              f"{v['nominal_empirical']*100:5.2f}% , delivers "
              f"{v['realised_empirical']*100:5.2f}%  ({v['empirical_ratio']:5.1f}x)"
              f"   threshold z > {v['threshold_median']:.2f}")

    print(f"\ndetection rate against an injected shift in {PRIMARY} "
          f"(robust SD units), at a 1% calibrated false-alarm rate")
    for sname in ("single", "scan"):
        print(f"  [{sname}]")
        print(f"{'observation':>18s} " + "".join(f"{e:8.2f}" for e in EFFECTS))
        row = [out["power"][f"{sname}|{e}"] for e in EFFECTS]
        print(f"{'1 record (0.4 min)':>18s} " + "".join(f"{v:8.2f}" for v in row))
        for k in POOL:
            row = [out["pooling"].get(f"{sname}|{k}|{e}") for e in EFFECTS]
            if any(v is None for v in row):
                continue
            mins = k * out["config"]["median_duration_s"] / 60.0
            print(f"{f'{k} records ({mins:.0f} min)':>18s} "
                  + "".join(f"{v:8.2f}" for v in row))


if __name__ == "__main__":
    main()
