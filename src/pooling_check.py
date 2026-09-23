"""Read the pooled-power table at zero effect: what false-alarm rate the
"calibrated 1 %" detector actually delivers once k recordings are averaged.

Uses only results/exp2_detector.json and results/exp1_variability.json (no
recordings needed).  Writes results/pooling_check.json and prints the tables
that went into the README section "Where the pooled power table gets it
wrong".
"""

import json
import math

ALPHA = 0.01


def main():
    det = json.load(open("results/exp2_detector.json"))
    var = json.load(open("results/exp1_variability.json"))
    pool = det["config"]["pool"]
    effects = det["config"]["effects_sd"]
    P = det["pooling"]

    out = {"alpha": ALPHA, "pool": pool, "false_alarm_vs_k": {},
           "net_power_1sd": {}, "root_k_prediction": {}, "tails_2_3_4_sigma": {}}

    # --- false-alarm rate of the calibrated detector as a function of k
    for s in ("single", "scan"):
        rows = []
        for k in pool:
            fa = P[f"{s}|{k}|0.0"]
            rows.append({"k": k, "minutes": round(k * det["config"]["median_duration_s"] / 60, 1),
                         "false_alarm": fa, "ratio_to_nominal": fa / ALPHA})
        out["false_alarm_vs_k"][s] = rows

    # --- power at 1 SD net of the false-alarm rate at the same k
    for s in ("single", "scan"):
        out["net_power_1sd"][s] = {
            str(k): {"power": P[f"{s}|{k}|1.0"], "false_alarm": P[f"{s}|{k}|0.0"],
                     "net": P[f"{s}|{k}|1.0"] - P[f"{s}|{k}|0.0"]} for k in pool}

    # --- the sqrt(k) mechanism.  A mean of k records shrinks the within-fish
    # spread by sqrt(k); the offset between the training fish and the held-out
    # fish does not shrink.  With ICC rho for the primary metric the share of
    # between-fish variance in a mean of k is rho_k = k rho / (1 + (k-1) rho).
    rho = var["icc"][det["false_alarms"]["primary_metric"]]["icc"]
    out["root_k_prediction"] = {
        "primary_metric": det["false_alarms"]["primary_metric"],
        "icc": rho,
        "between_share_of_mean_of_k": {str(k): k * rho / (1 + (k - 1) * rho) for k in pool}}

    # --- bout-level tails at 2, 3 and 4 sigma
    for q, t in var["tails_bout_level"].items():
        out["tails_2_3_4_sigma"][q] = {
            "excess_kurtosis": t["excess_kurtosis"],
            "ratio_2sigma": t["2sigma"]["ratio"],
            "ratio_3sigma": t["3sigma"]["ratio"],
            "ratio_4sigma": t["4sigma"]["ratio"],
            "realised_4sigma_pct": 100 * t["4sigma"]["realised"],
            "one_in_4sigma": (round(1 / t["4sigma"]["realised"]) if t["4sigma"]["realised"] > 0 else None)}

    with open("results/pooling_check.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"calibrated at alpha = {ALPHA}; false-alarm rate on held-out fish vs k\n")
    print("k   min   single   ratio   scan    ratio")
    for a, b in zip(out["false_alarm_vs_k"]["single"], out["false_alarm_vs_k"]["scan"]):
        print(f"{a['k']:>2} {a['minutes']:>5} {100*a['false_alarm']:6.2f}%  {a['ratio_to_nominal']:4.1f}x "
              f"{100*b['false_alarm']:6.2f}%  {b['ratio_to_nominal']:4.1f}x")
    print(f"\nICC({out['root_k_prediction']['primary_metric']}) = {rho:.3f}; between-fish share of a mean of k:")
    print("  " + "  ".join(f"k={k}: {v:.2f}" for k, v in out["root_k_prediction"]["between_share_of_mean_of_k"].items()))
    print("\npower at 1 SD, net of false alarms")
    for s in ("single", "scan"):
        print(f"  {s:6s} " + "  ".join(f"k={k}: {v['power']:.2f}-{v['false_alarm']:.2f}={v['net']:.2f}"
                                       for k, v in out["net_power_1sd"][s].items()))
    print("\nbout-level tails: exceedance / Gaussian promise")
    print("quantity               kurt   2s     3s     4s    1-in at 4s")
    for q, t in out["tails_2_3_4_sigma"].items():
        print(f"{q:22s} {t['excess_kurtosis']:5.0f} {t['ratio_2sigma']:5.2f} {t['ratio_3sigma']:5.1f} "
              f"{t['ratio_4sigma']:6.0f}  {t['one_in_4sigma']}")


if __name__ == "__main__":
    main()
