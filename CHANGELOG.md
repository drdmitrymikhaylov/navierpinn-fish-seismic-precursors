# Changelog

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
