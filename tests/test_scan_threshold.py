"""Pin the section "Who sets the scan's threshold" (added 2026-10-19) to
results/scan_threshold_check.json and results/exp2_detector.json.  Needs no
recordings.  The table rows on the page are built by the formatters below
and pasted verbatim, so a rerun that moves a number fails here.
Run:  /usr/local/bin/python3 -m pytest tests/test_scan_threshold.py
"""

import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NAMES = {
    "bout_rate": "bout rate",
    "ibi_median": "inter-bout interval median",
    "ibi_cv": "inter-bout interval CV",
    "displacement_median": "displacement median",
    "distance_per_s": "distance per second",
    "turn_dispersion": "turn dispersion",
    "turn_abs_median": "median absolute turn",
    "large_turn_fraction": "large-turn fraction",
    "turn_persistence": "turn persistence",
    "centre_distance_median": "centre distance (thigmotaxis)",
    "area_explored": "area explored",
}
EFFECTS = ["0.0", "0.25", "0.5", "1.0", "1.5", "2.0", "3.0"]


def _load(name):
    with open(os.path.join(HERE, "results", name)) as fh:
        return json.load(fh)


def _readme():
    with open(os.path.join(HERE, "README.md"), encoding="utf-8") as fh:
        return " ".join(fh.read().split())


def tail_rows(d):
    t = d["per_metric_training_tail"]
    order = sorted(t, key=lambda n: -t[n]["q99_median"])
    return [f"| {NAMES[n]} | {t[n]['q99_median']:.2f} | {t[n]['q999_median']:.2f} "
            f"| {100 * t[n]['share_of_heldout_alarms']:.0f} % |" for n in order]


def variant_rows(d):
    loo, ref, cal = d["leave_one_metric_out"], d["reference"], d["per_metric_calibrated_scan"]
    stored = _load("exp2_detector.json")
    rows = [f"| all eleven | {d['scan_threshold_median']:.2f} | "
            f"{100 * stored['false_alarms']['scan']['realised_empirical']:.2f} % | "
            f"{ref['scan_k1_3sd']:.3f} | {ref['scan_k16_1sd']:.3f} |"]
    named = ["turn_dispersion", "turn_abs_median", "ibi_cv"]
    for n in named:
        v = loo[n]
        rows.append(f"| without {NAMES[n]} | {v['threshold_median']:.2f} | "
                    f"{100 * v['fa_delivered']:.2f} % | {v['power_k1_3sd']:.3f} | "
                    f"{v['power_k16_1sd']:.3f} |")
    rest = [loo[n] for n in loo if n not in named and n != d["primary_metric"]]
    assert len(rest) == 7
    lo = lambda k: min(v[k] for v in rest)   # noqa: E731
    hi = lambda k: max(v[k] for v in rest)   # noqa: E731
    rows.append(f"| without any one of the other seven | "
                f"{lo('threshold_median'):.2f}–{hi('threshold_median'):.2f} | "
                f"{100 * lo('fa_delivered'):.2f}–{100 * hi('fa_delivered'):.2f} % | "
                f"{lo('power_k1_3sd'):.3f}–{hi('power_k1_3sd'):.3f} | "
                f"{lo('power_k16_1sd'):.3f}–{hi('power_k16_1sd'):.3f} |")
    rows.append(f"| each metric in its own tail (q = {cal['q_median']:.4f}) | "
                f"{cal['primary_threshold_median']:.2f} on bout rate | "
                f"{100 * cal['fa_delivered_k1']:.2f} % | {cal['power_k1_3sd']:.3f} | "
                f"{cal['power_k16_1sd']:.3f} at {100 * cal['fa_delivered_k16']:.1f} % |")
    rows.append(f"| single metric, bout rate | {d['single_threshold_median']:.2f} | "
                f"{100 * stored['false_alarms']['single']['realised_empirical']:.2f} % | "
                f"{ref['single_k1_3sd']:.3f} | {ref['single_k16_1sd']:.3f} |")
    return rows


def mde_rows(d):
    med = _load("exp2_detector.json")["config"]["median_duration_s"]
    f = lambda x: "> 3" if x is None else f"{x:.2f}"   # noqa: E731
    return [f"| {k} | {int(k) * med / 60:.1f} | {f(d['mde_80pct_sd']['single'][k])} "
            f"| {f(d['mde_80pct_sd']['scan'][k])} |" for k in ["1", "2", "4", "8", "16"]]


def test_reproduces_exp2():
    d = _load("scan_threshold_check.json")
    assert d["reproduces_exp2_pooling_to"] == 0.0
    assert d["n_splits"] == 200 == _load("exp2_detector.json")["false_alarms"]["n_splits"]


def test_thresholds():
    d = _load("scan_threshold_check.json")
    s = _load("exp2_detector.json")["false_alarms"]
    assert round(d["scan_threshold_median"], 2) == 11.97
    assert round(s["scan"]["threshold_median"], 2) == 11.97
    assert round(d["single_threshold_median"], 2) == 3.15
    assert round(d["scan_threshold_median"] / d["single_threshold_median"], 1) == 3.8
    r = _readme()
    assert "robust z of **11.97**" in r and "**3.15** for the single metric" in r
    assert "3.8 times higher" in r


def test_who_owns_the_threshold():
    d = _load("scan_threshold_check.json")
    t = d["per_metric_training_tail"]
    share = {n: round(100 * t[n]["share_of_heldout_alarms"]) for n in t}
    assert share["turn_dispersion"] == 57 and share["turn_abs_median"] == 28
    assert share["ibi_cv"] == 9 and share["ibi_median"] == 6
    zero = [n for n in t if t[n]["share_of_heldout_alarms"] == 0]
    assert len(zero) == 7 and "bout_rate" in zero
    assert abs(sum(t[n]["share_of_heldout_alarms"] for n in t) - 1) < 1e-12
    # the scan threshold sits above every metric's own 99.9th percentile but one
    above = [n for n in t if t[n]["q999_median"] > d["scan_threshold_median"]]
    assert above == ["turn_dispersion"]
    q999_zero = [t[n]["q999_median"] for n in zero]
    assert round(min(q999_zero), 1) == 2.9 and round(max(q999_zero), 1) == 5.8
    r = _readme()
    for row in tail_rows(d):
        assert row in r, row
    assert "arg-max in **57 %**" in r and "median absolute turn in 28 %" in r
    assert "seven of the eleven metrics" in r


def test_scan_is_exactly_blind_at_one_recording():
    d = _load("scan_threshold_check.json")
    p = _load("exp2_detector.json")["pooling"]
    v = [d["scan_power_equals_fa_at_k1"][e] for e in EFFECTS]
    assert max(v) - min(v) < 5e-5 and round(100 * v[0], 2) == 1.50
    assert v[0] == p["scan|1|0.0"]
    assert round(100 * d["scan_power_k2"]["3.0"], 1) == 4.1
    assert round(100 * p["scan|4|3.0"]) == 63
    r = _readme()
    assert "**1.50 % at 0, 0.25, 0.5, 1, 1.5, 2 and 3 SD alike**" in r
    assert "At two recordings a 3 SD shift lifts it to 4.1 %; at four, 63 %" in r


def test_leave_one_out_and_calibrated_scan():
    d = _load("scan_threshold_check.json")
    loo, cal = d["leave_one_metric_out"], d["per_metric_calibrated_scan"]
    assert round(loo["turn_dispersion"]["threshold_median"], 2) == 9.73
    gains = [loo[n]["power_k16_1sd"] - d["reference"]["scan_k16_1sd"]
             for n in loo if n != d["primary_metric"]]
    assert round(max(gains), 2) == 0.04 and min(gains) > 0
    assert round(loo["turn_dispersion"]["power_k16_1sd"]
                 - d["reference"]["scan_k16_1sd"], 2) == 0.03
    assert round(cal["q_median"], 4) == 0.9998
    assert round(100 * cal["fa_delivered_k1"], 1) == 3.8
    assert round(100 * cal["fa_delivered_k16"], 1) == 10.9
    assert round(cal["power_k1_3sd"], 2) == 0.23 and round(cal["power_k16_1sd"], 2) == 0.69
    assert round(cal["primary_threshold_median"], 2) == 3.97
    r = _readme()
    for row in variant_rows(d):
        assert row in r, row
    assert "**0.9998**" in r and "delivers **3.8 %** false alarms" in r


def test_minimum_detectable_effect():
    d = _load("scan_threshold_check.json")
    m = d["mde_80pct_sd"]
    assert m["single"]["1"] is None and m["scan"]["4"] is None
    assert round(m["single"]["16"], 2) == 0.97 and round(m["scan"]["16"], 2) == 1.35
    assert round(m["single"]["4"], 2) == 1.91 and round(m["scan"]["8"], 2) == 2.04
    r = _readme()
    for row in mde_rows(d):
        assert row in r, row
