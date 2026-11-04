# Binding at the rank thresholds (L = 1, n = 4)

`scripts/threshold_binding.py`

| target | d | N_V* | log2 q | log2 beta | binding (min over models) | margin | S5 feasible | log2 budget | log2 max sigma |
|---|---|---|---|---|---|---|---|---|---|
| 128 classical | 6 | 180 | 38.544 | 38.544 | 130.7 | 2.7 | True | -171.98 | 20.62 |
| 128 classical | 7 | 3021514 | 45.250 | 45.250 | 130.1 | 2.1 | True | -170.68 | 27.58 |
| 128 classical+quantum | 7 | 354255 | 43.745 | 43.745 | 130.4 | 2.4 | True | -171.23 | 26.03 |
| 128 core-SVP | 8 | 93635 | 43.026 | 43.026 | 132.6 | 4.6 | True | -171.74 | 25.08 |
| 128 core-SVP classical+quantum | 8 | 900 | 39.985 | 39.985 | 133.3 | 5.3 | True | -171.98 | 21.86 |
| 128 core-SVP classical+quantum | 9 | 1852773 | 45.162 | 45.162 | 132.8 | 4.8 | True | -170.50 | 27.25 |
| 192 classical | 9 | 420 | 39.585 | 39.585 | 200.8 | 8.8 | True | -171.98 | 21.38 |
| 192 classical | 10 | 400666 | 44.226 | 44.226 | 199.9 | 7.9 | True | -171.15 | 26.15 |
| 192 classical+quantum | 10 | 9400 | 41.913 | 41.913 | 200.3 | 8.3 | True | -171.95 | 23.64 |
| 192 classical+quantum | 11 | 6846195 | 46.315 | 46.315 | 199.4 | 7.4 | True | -170.57 | 28.21 |
| hiding 192 / binding 128, classical | 9 | 420 | 39.585 | 39.585 | 200.8 | 72.8 | True | -171.98 | 21.38 |
| hiding 192 / binding 128, classical | 10 | 400666 | 44.226 | 44.226 | 199.9 | 71.9 | True | -171.15 | 26.15 |
