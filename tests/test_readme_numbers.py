"""Pin the numbers quoted in README.md to results/*.json, so the page cannot
drift from the files.  Needs no recordings: everything is read from results/.
Run:  /usr/local/bin/python3 -m pytest tests/test_readme_numbers.py
"""

import json
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name):
    with open(os.path.join(HERE, "results", name)) as fh:
        return json.load(fh)


def _readme():
    with open(os.path.join(HERE, "README.md"), encoding="utf-8") as fh:
        return fh.read()


def test_dataset_counts():
    v = _load("exp1_variability.json")["config"]
    r = _readme()
    assert v["n_recordings"] == 504 and "504 recordings survive" in r
    assert v["n_fish"] == 39 and "39 larval zebrafish" in r
    assert round(v["median_duration_s"]) == 26
    assert round(v["total_hours"], 1) == 6.4 and "6.4 hours" in r


def test_bout_level_tails_2_3_4_sigma():
    t = _load("exp1_variability.json")["tails_bout_level"]
    r = _readme()
    # kurtosis and 3-sigma exceedance as printed
    assert round(t["turn_angle_deg"]["excess_kurtosis"], 1) == 47.6
    assert round(t["inter_bout_interval_s"]["excess_kurtosis"], 1) == 414.9
    assert round(100 * t["turn_angle_deg"]["3sigma"]["realised"], 2) == 1.86
    assert round(100 * t["inter_bout_interval_s"]["3sigma"]["realised"], 2) == 1.09
    # the 2-sigma / 4-sigma columns added 2026-09-23
    assert round(t["turn_angle_deg"]["2sigma"]["ratio"], 1) == 1.1
    assert round(t["inter_bout_interval_s"]["2sigma"]["ratio"], 1) == 0.4
    assert round(t["turn_angle_deg"]["4sigma"]["ratio"]) == 133
    assert round(t["displacement"]["4sigma"]["ratio"]) == 49
    assert round(t["inter_bout_interval_s"]["4sigma"]["ratio"]) == 102
    assert round(1 / t["turn_angle_deg"]["4sigma"]["realised"]) == 119
    assert "**133×** (1 in 119)" in r
    # 3 -> 4 sigma: promised drops 43x, delivered 2.2x
    promised = t["turn_angle_deg"]["3sigma"]["nominal"] / t["turn_angle_deg"]["4sigma"]["nominal"]
    delivered = t["turn_angle_deg"]["3sigma"]["realised"] / t["turn_angle_deg"]["4sigma"]["realised"]
    assert round(promised) == 43 and round(delivered, 1) == 2.2
    assert round(1 / t["turn_angle_deg"]["4sigma"]["nominal"], -2) == 15800


def test_icc_table_and_ceiling():
    icc = _load("exp1_variability.json")["icc"]
    assert max(v["icc"] for v in icc.values()) < 0.285
    assert round(icc["displacement_median"]["icc"], 2) == 0.28
    assert round(icc["turn_dispersion"]["icc"], 2) == 0.13
    assert round(icc["centre_distance_median"]["icc"], 2) == 0.08
    assert round(icc["turn_persistence"]["icc"], 2) == 0.02
    assert round(icc["bout_rate"]["icc"], 2) == 0.22


def test_frozen_detector_table():
    f = _load("exp2_detector.json")["false_alarms"]
    assert f["n_splits"] == 200 and f["primary_metric"] == "bout_rate"
    assert round(100 * f["single"]["realised_gaussian"], 2) == 1.42
    assert round(100 * f["single"]["realised_empirical"], 2) == 1.53
    assert round(100 * f["scan"]["realised_gaussian"], 2) == 35.64
    assert round(100 * f["scan"]["realised_empirical"], 2) == 1.50
    assert round(f["scan"]["gaussian_ratio"], 1) == 12.2


def test_power_tables():
    d = _load("exp2_detector.json")
    P, S = d["pooling"], d["power"]
    one = [round(S[f"single|{e}"], 2) for e in (0.5, 1.0, 1.5, 2.0, 3.0)]
    assert one == [0.02, 0.05, 0.09, 0.17, 0.47]
    assert [round(P[f"single|16|{e}"], 2) for e in (0.5, 1.0, 1.5)] == [0.29, 0.83, 0.99]
    assert [round(P[f"scan|16|{e}"], 2) for e in (1.0, 1.5, 2.0)] == [0.51, 0.93, 1.00]
    assert [round(P[f"scan|8|{e}"], 2) for e in (1.0, 1.5, 2.0)] == [0.12, 0.40, 0.79]


def test_pooled_false_alarm_rate_rises_with_k():
    """The 2026-09-23 finding: the calibrated 1 % detector delivers 4.0 %
    (single) and 6.7 % (scan) once 16 records are pooled."""
    d = _load("exp2_detector.json")
    P = d["pooling"]
    c = _load("pooling_check.json")
    r = _readme()
    single = {row["k"]: row["false_alarm"] for row in c["false_alarm_vs_k"]["single"]}
    scan = {row["k"]: row["false_alarm"] for row in c["false_alarm_vs_k"]["scan"]}
    for k in d["config"]["pool"]:
        assert single[k] == P[f"single|{k}|0.0"]
        assert scan[k] == P[f"scan|{k}|0.0"]
    fa_s = [round(100 * single[k], 2) for k in (1, 2, 4, 8, 16)]
    fa_c = [round(100 * scan[k], 2) for k in (1, 2, 4, 8, 16)]
    assert fa_s == [1.58, 1.63, 1.91, 2.38, 4.02]
    assert fa_c == [1.50, 1.85, 2.24, 3.54, 6.66]
    # monotone in k for both statistics
    assert all(a < b for a, b in zip(fa_s, fa_s[1:]))
    assert all(a < b for a, b in zip(fa_c, fa_c[1:]))
    assert "| 16 | 7.0 | **4.02%** | 4.0× | **6.66%** | 6.7× |" in r
    # net power at 1 SD, k = 16
    net = c["net_power_1sd"]
    assert round(net["single"]["16"]["net"], 2) == 0.79
    assert round(net["scan"]["16"]["net"], 2) == 0.44
    # the ICC arithmetic quoted in the text
    share = c["root_k_prediction"]["between_share_of_mean_of_k"]
    rho = c["root_k_prediction"]["icc"]
    assert round(rho, 2) == 0.22
    for k, v in share.items():
        k = int(k)
        assert abs(v - k * rho / (1 + (k - 1) * rho)) < 1e-12
    assert round(share["1"], 2) == 0.22 and round(share["4"], 2) == 0.53
    assert round(share["16"], 2) == 0.82
    assert re.search(r"\*\*0\.82 at\s+sixteen\*\*", r)


def test_space_use_numbers():
    e = _load("exp3_fp_pinn.json")
    assert e["n_transitions"] == 15450
    assert e["relaxation"]["tau_fd_s"] == 155.0 == e["relaxation"]["tau_pinn_s"]
    assert round(100 * e["relaxation"]["pinn_vs_fd_max_rel_err"], 1) == 1.4
    assert round(e["stationary"]["kl_obs_vs_model"], 1) == 4.4
    assert round(e["stationary"]["median_r_model"], 2) == 0.96
    assert round(e["stationary"]["median_r_observed"], 2) == 0.55  # was printed 0.56 until 2026-09-23
    assert "observed\nmedian is 0.55" in _readme() or "observed median is 0.55" in _readme()
    assert round(e["recordings"]["frac_last_bout_at_wall"], 2) == 0.89
    assert round(e["stress_response"]["effect_size_sd_at_26s"]["1.5"], 2) == 0.66
    assert round(e["stress_response"]["effect_size_sd_at_26s"]["2.0"], 2) == 0.85


def test_pooling_check_json_matches_readme_minutes():
    c = _load("pooling_check.json")
    med = _load("exp2_detector.json")["config"]["median_duration_s"]
    for row in c["false_alarm_vs_k"]["single"]:
        assert row["minutes"] == round(row["k"] * med / 60, 1)
    assert [row["minutes"] for row in c["false_alarm_vs_k"]["single"]] == [0.4, 0.9, 1.7, 3.5, 7.0]
