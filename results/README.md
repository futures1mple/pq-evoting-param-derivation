# results

| Path | Written by | Content |
|---|---|---|
| `estimates/lwe.jsonl` | `estimate_grid.py`, `evolve_chain.py estimate` | M-LWE estimates at q_b, d = 6..12 and b = 36..50 (420 entries, `estimate_grid.py`), plus d = 5 and b = 35 (44 entries, `evolve_chain.py`); four cost models |
| `estimates/sis.jsonl` | `estimate_grid.py` | M-SIS estimates of T-EVOLVE at q_b (1856 entries) |
| `estimates/lwe_fine.jsonl`, `crossings.json` | `ceilings.py` | M-LWE at +-0.05 bit around each interpolated crossing of 128 and 192 bits, and the midpoint refinements (`--refine`, `--refine-dual-free`) that bring every ceiling to 0.1 bit |
| `estimates/sis_fine.jsonl` | `derive.py estimate` | M-SIS at the selected moduli |
| `estimates/lwe_sel.jsonl`, `estimates/lwe_L2.jsonl` | `derive.py estimate` | M-LWE at the selected moduli, L = 1 and L = 2 |
| `estimates/lwe_naive.jsonl` | `naive_scaling.py` | M-LWE at the moduli of the simpler strategies |
| `estimates/sis_thr.jsonl` | `threshold_binding.py` | M-SIS at the rank thresholds |
| `estimates/bounds.jsonl` | `estimate_grid.py` | bounds of the chain at the grid configurations |
| `estimates/{SUMMARY,PROGRESS,machine}.txt`, `run.log`, `run_info.jsonl` | `estimate_grid.py` | logs of the main run (server, 7 processes) |
| `estimates_evolve/sis.jsonl` | `evolve_chain.py estimate` | M-SIS estimates of the EVOLVE module (1568 entries) |
| `derivation.{json,md}` | `derive.py select` | selected sets for every target, hiding ceilings, thresholds, ceilings without dual attacks |
| `threshold_binding.{json,md}` | `threshold_binding.py` | binding and feasibility at the thresholds |
| `naive.{json,md}` | `naive_scaling.py` | simpler strategies at 6.48 million voters |
| `sweep_electorate.json`, `sweep_d8.json` | `sweep_electorate.py`, `figures.py` | bounds on a dense electorate grid (rank 7; rank 8 for the figure) |
| `paper_numbers.txt` | `paper_numbers.py` | table rows and derived numbers of the paper (Tables 2-5; effect of each input, authority-side costs and T-EVOLVE vs EVOLVE bounds in Section 6 and 7) |
| `evolve_selection.{json,md}` | `evolve_chain.py select` | superseded: EVOLVE selection restricted to the moduli q_b (the paper uses `derivation.json`) |
| `selection.{json,md}` | `select_params.py` | superseded: T-EVOLVE selection restricted to the moduli q_b |
| `figures/` | `figures.py` | Figure 1 of the paper and the rank plot |

Every estimate line contains the instance, the cost model, the cost of every attack, the block sizes,
the minimum and the estimator commit.
