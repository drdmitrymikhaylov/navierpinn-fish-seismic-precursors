"""Who owns the scan's threshold, and what it costs.

The scan statistic of exp2_detector.py is the largest absolute robust z over
eleven metrics, thresholded at the 99th percentile of the training statistic.
results/exp2_detector.json stores that threshold's median as 11.97 -- four
times the single-metric threshold of 3.15 -- and the stored power of the scan
at one recording equals its false-alarm rate, to four decimals, at every
injected shift up to 3 SD.  This script asks why, on the same 200 splits of
fish as exp2 (same seeds, same training fraction, same pooled draws):

  1. per metric, the 99th and 99.9th percentile of |z| on the training
     recordings, and which metric is the arg-max of each held-out scan
     exceedance -- who sets the threshold and who fires;
  2. leave-one-metric-out: the scan over the other ten, its threshold,
     delivered false alarms, and power at (k = 1, 3 SD) and (k = 16, 1 SD);
  3. a per-metric calibrated scan: one quantile q for every metric, chosen
     so that the union of the eleven exceedances is 1 % of the training
     records (so each metric is tested in its own 1-in-1100 tail);
  4. the minimum detectable effect at 80 % power in the stored power tables,
     by linear interpolation between grid points (stored JSON only).

The exp2 pooled sweep is recomputed on the way and must agree with the stored
numbers; the script asserts it.  Writes results/scan_threshold_check.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from behaviour import load, state_matrix, valid          # noqa: E402
from exp2_detector import (ALPHA, EFFECTS, POOL, N_SPLITS, PRIMARY,  # noqa
                           robust_scale, stat_scan, stat_single, zscores)

CHECK_SPLITS = 200


def union_quantile(Ztr, alpha):
    """Largest common per-metric quantile q whose union exceedance on the
    training records is <= alpha; returns q and the per-metric thresholds."""
    n, m = Ztr.shape
    # candidate q on a fine grid between 1 - alpha and 1 - alpha / (m * 4)
    best = None
    for q in np.linspace(1 - alpha, 1 - alpha / (4 * m), 400):
        thr = np.quantile(Ztr, q, axis=0)
        rate = float(np.mean((Ztr > thr).any(axis=1)))
        if rate <= alpha:
            best = (float(q), thr, rate)
            break
    if best is None:
        q = 1 - alpha / (4 * m)
        thr = np.quantile(Ztr, q, axis=0)
        best = (float(q), thr, float(np.mean((Ztr > thr).any(axis=1))))
    return best


def mde(power_by_effect, target=0.8):
    """Smallest effect (SD) reaching `target` power, linear between grid
    points; None if the largest grid effect does not reach it."""
    es = sorted(power_by_effect)
    for a, b in zip(es, es[1:]):
        pa, pb = power_by_effect[a], power_by_effect[b]
        if pa < target <= pb:
            return float(a + (target - pa) / (pb - pa) * (b - a))
    return None


def main():
    trials = [t for t in load() if valid(t)]
    X, fish, names = state_matrix(trials)
    keep = np.isfinite(X).all(axis=1)
    X, fish = X[keep], fish[keep]
    fish_ids = np.unique(fish)
    m = len(names)
    jp = names.index(PRIMARY)
    E1, E16 = 3.0, 1.0            # the two shifts read for every variant

    q99 = {n: [] for n in names}
    q999 = {n: [] for n in names}
    owner_test = {n: 0 for n in names}      # arg-max metric of held-out alarms
    owner_train = {n: 0 for n in names}     # arg-max metric of training top 1 %
    scan_thr = []
    rep = {f"{s}|{k}|{e}": [] for s in ("scan", "single")
           for k in POOL for e in EFFECTS}
    loo = {n: {"thr": [], "fa": [], "p1_3sd": [], "p16_1sd": []} for n in names}
    cal = {"q": [], "fa": [], "p1_3sd": [], "p16_1sd": [],
           "fa16": [], "thr_primary": []}
    n_used = 0

    for split in range(CHECK_SPLITS):
        rng = np.random.default_rng(split)
        perm = rng.permutation(fish_ids)
        tr_fish = set(perm[: int(0.6 * len(perm))].tolist())
        tr = np.array([f in tr_fish for f in fish])
        if tr.sum() < 20 or (~tr).sum() < 20:
            continue
        n_used += 1
        med, mad = robust_scale(X[tr])
        Xte = X[~tr]
        Ztr = np.abs(zscores(X[tr], med, mad))
        Zte = np.abs(zscores(Xte, med, mad))

        # 1. who owns the threshold
        t_scan = float(np.quantile(Ztr.max(axis=1), 1 - ALPHA))
        scan_thr.append(t_scan)
        for j, n in enumerate(names):
            q99[n].append(float(np.quantile(Ztr[:, j], 0.99)))
            q999[n].append(float(np.quantile(Ztr[:, j], 0.999)))
        for i in np.where(Zte.max(axis=1) > t_scan)[0]:
            owner_test[names[int(np.argmax(Zte[i]))]] += 1
        for i in np.where(Ztr.max(axis=1) > t_scan)[0]:
            owner_train[names[int(np.argmax(Ztr[i]))]] += 1

        # 2./3. variants at k = 1
        t_single = float(np.quantile(Ztr[:, jp], 1 - ALPHA))
        q, thr_cal, _ = union_quantile(Ztr, ALPHA)
        cal["q"].append(q)
        cal["thr_primary"].append(float(thr_cal[jp]))
        cal["fa"].append(float(np.mean((Zte > thr_cal).any(axis=1))))
        Z = Xte.copy()
        Z[:, jp] += E1 * mad[jp]
        Zs = np.abs(zscores(Z, med, mad))
        cal["p1_3sd"].append(float(np.mean((Zs > thr_cal).any(axis=1))))
        for j, n in enumerate(names):
            cols = [c for c in range(m) if c != j]
            t_j = float(np.quantile(Ztr[:, cols].max(axis=1), 1 - ALPHA))
            loo[n]["thr"].append(t_j)
            loo[n]["fa"].append(float(np.mean(Zte[:, cols].max(axis=1) > t_j)))
            loo[n]["p1_3sd"].append(
                float(np.mean(Zs[:, cols].max(axis=1) > t_j))
                if j != jp else float("nan"))

        # pooling, drawn in the same order as exp2 so the sets are identical
        for k in POOL:
            n_draw = 400
            gt = rng.integers(0, tr.sum(), size=(n_draw, k))
            base_tr = X[tr][gt].mean(axis=1)
            gte = rng.integers(0, len(Xte), size=(n_draw, k))
            base = Xte[gte].mean(axis=1)
            med_k, mad_k = robust_scale(base_tr)
            Btr = np.abs(zscores(base_tr, med_k, mad_k))
            t_k = {"scan": float(np.quantile(Btr.max(axis=1), 1 - ALPHA)),
                   "single": float(np.quantile(Btr[:, jp], 1 - ALPHA))}
            for e in EFFECTS:
                Z = base.copy()
                Z[:, jp] += e * mad[jp]
                rep[f"scan|{k}|{e}"].append(
                    float(np.mean(stat_scan(Z, med_k, mad_k) > t_k["scan"])))
                rep[f"single|{k}|{e}"].append(
                    float(np.mean(stat_single(Z, med_k, mad_k, jp)
                                  > t_k["single"])))
            if k == 16:
                Bte = np.abs(zscores(base, med_k, mad_k))
                Z = base.copy()
                Z[:, jp] += E16 * mad[jp]
                Bs = np.abs(zscores(Z, med_k, mad_k))
                qk, thr_k, _ = union_quantile(Btr, ALPHA)
                cal["fa16"].append(float(np.mean((Bte > thr_k).any(axis=1))))
                cal["p16_1sd"].append(float(np.mean((Bs > thr_k).any(axis=1))))
                for j, n in enumerate(names):
                    cols = [c for c in range(m) if c != j]
                    t_j = float(np.quantile(Btr[:, cols].max(axis=1), 1 - ALPHA))
                    loo[n]["p16_1sd"].append(
                        float(np.mean(Bs[:, cols].max(axis=1) > t_j))
                        if j != jp else float("nan"))

    stored = json.load(open("results/exp2_detector.json"))
    rep_mean = {k: float(np.mean(v)) for k, v in rep.items() if v}
    worst = max(abs(rep_mean[k] - stored["pooling"][k]) for k in rep_mean)
    assert worst < 1e-9, f"pooled sweep does not reproduce: {worst}"
    assert n_used == stored["false_alarms"]["n_splits"]

    n_alarms = sum(owner_test.values())
    out = {
        "n_splits": n_used, "alpha": ALPHA, "primary_metric": PRIMARY,
        "reproduces_exp2_pooling_to": worst,
        "scan_threshold_median": float(np.median(scan_thr)),
        "single_threshold_median": stored["false_alarms"]["single"]["threshold_median"],
        "per_metric_training_tail": {
            n: {"q99_median": float(np.median(q99[n])),
                "q999_median": float(np.median(q999[n])),
                "share_of_heldout_alarms": owner_test[n] / n_alarms,
                "share_of_training_top1pct": owner_train[n] / sum(owner_train.values())}
            for n in names},
        "heldout_alarms_total": n_alarms,
        "scan_power_equals_fa_at_k1": {
            e: stored["pooling"][f"scan|1|{e}"] for e in EFFECTS},
        "scan_power_k2": {e: stored["pooling"][f"scan|2|{e}"] for e in EFFECTS},
        "leave_one_metric_out": {
            n: {"threshold_median": float(np.median(loo[n]["thr"])),
                "fa_delivered": float(np.mean(loo[n]["fa"])),
                "power_k1_3sd": (None if n == PRIMARY
                                 else float(np.mean(loo[n]["p1_3sd"]))),
                "power_k16_1sd": (None if n == PRIMARY
                                  else float(np.mean(loo[n]["p16_1sd"])))}
            for n in names},
        "per_metric_calibrated_scan": {
            "q_median": float(np.median(cal["q"])),
            "primary_threshold_median": float(np.median(cal["thr_primary"])),
            "fa_delivered_k1": float(np.mean(cal["fa"])),
            "power_k1_3sd": float(np.mean(cal["p1_3sd"])),
            "fa_delivered_k16": float(np.mean(cal["fa16"])),
            "power_k16_1sd": float(np.mean(cal["p16_1sd"]))},
        "reference": {
            "scan_k1_3sd": stored["pooling"]["scan|1|3.0"],
            "scan_k16_1sd": stored["pooling"]["scan|16|1.0"],
            "single_k1_3sd": stored["pooling"]["single|1|3.0"],
            "single_k16_1sd": stored["pooling"]["single|16|1.0"]},
        "mde_80pct_sd": {
            s: {k: mde({e: stored["pooling"][f"{s}|{k}|{e}"] for e in EFFECTS})
                for k in POOL} for s in ("single", "scan")},
    }
    os.makedirs("results", exist_ok=True)
    with open("results/scan_threshold_check.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"splits {n_used}; exp2 pooled sweep reproduced to {worst:.1e}")
    print(f"scan threshold (median) {out['scan_threshold_median']:.2f}, "
          f"single {out['single_threshold_median']:.2f}")
    print(f"{'metric':24s} {'q99':>6s} {'q99.9':>7s} {'alarms':>7s} {'top1%':>6s}")
    for n in names:
        t = out["per_metric_training_tail"][n]
        print(f"{n:24s} {t['q99_median']:6.2f} {t['q999_median']:7.2f} "
              f"{100*t['share_of_heldout_alarms']:6.1f}% "
              f"{100*t['share_of_training_top1pct']:5.1f}%")
    print("\nleave one metric out      thr     FA   p(k1,3SD)  p(k16,1SD)")
    for n in names:
        t = out["leave_one_metric_out"][n]
        p1 = "    (primary)" if n == PRIMARY else (
            f"{t['power_k1_3sd']:8.3f} {t['power_k16_1sd']:10.3f}")
        print(f"{n:24s} {t['threshold_median']:6.2f} "
              f"{100*t['fa_delivered']:5.2f}% {p1}")
    c = out["per_metric_calibrated_scan"]
    print(f"\nper-metric calibrated scan: q {c['q_median']:.4f}, primary thr "
          f"{c['primary_threshold_median']:.2f}, FA k1 {100*c['fa_delivered_k1']:.2f}%, "
          f"p(k1,3SD) {c['power_k1_3sd']:.3f}, FA k16 {100*c['fa_delivered_k16']:.2f}%, "
          f"p(k16,1SD) {c['power_k16_1sd']:.3f}")
    print("reference:", out["reference"])
    print("MDE at 80%:", out["mde_80pct_sd"])


if __name__ == "__main__":
    main()
