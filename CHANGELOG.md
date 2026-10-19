# Changelog

## 2026-10-19

- LICENSE now holds the standard MIT text and LICENSE-CONTENT the standard CC BY 4.0 text, so GitHub recognises the licence (it showed "Other" before). The terms are unchanged; the note on what each licence covers moved to `NOTICE.md`.
- The scan's calibrated threshold (robust z 11.97 vs 3.15 for the single
  metric) was never explained. Re-run of exp2's 200 splits with the arg-max
  kept: turn dispersion owns 57 % of the scan's alarms, median absolute turn
  28 %, seven of the eleven metrics (the primary among them) none. At one
  recording the scan's power equals its false-alarm rate (1.50 %) at every
  shift up to 3 SD — exactly blind, not "nearly". Leave-one-metric-out buys
  at most 0.04 of power; a per-metric calibrated scan needs each metric's
  sample maximum (q = 0.9998) and delivers 3.8 % / 10.9 % false alarms at
  1 / 16 recordings. Minimum detectable effect at 80 %: single 2.8 / 1.9 /
  1.4 / 0.97 SD at 2 / 4 / 8 / 16 recordings, scan 2.0 / 1.35 at 8 / 16.
  New section "Who sets the scan's threshold", finding 8, two sentences in
  "How long a look" qualified ("eleven heavy tails" → two).
  `src/scan_threshold_check.py`, `results/scan_threshold_check.json`,
  `tests/test_scan_threshold.py` (6 tests). Two tracked `._*` files removed.

## 2026-09-23

- The pooled power table was read at "the calibrated 1 % false-alarm rate";
  the zero-effect column of the same sweep says the pooled detector delivers
  1.6 → 4.0 % (single) and 1.5 → 6.7 % (scan) as k goes 1 → 16. New section
  "Where the pooled power table gets it wrong" with the table, the √k / ICC
  arithmetic (between-fish share 0.22 → 0.82 at k = 16) and the corrected
  reading of findings 5 and 6. Bout-level tail table gets 2σ and 4σ columns
  (turn angle: 133× the Gaussian promise at 4σ, one in 119). Observed median
  radius corrected 0.56 → 0.55; "pooled-power sweeps are held privately" was
  wrong (it is in `exp2_detector.py`). Added `src/pooling_check.py`,
  `results/pooling_check.json`, `tests/test_readme_numbers.py` (8 tests).
