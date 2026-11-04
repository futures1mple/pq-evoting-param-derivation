# pq-evoting-param-derivation

Parameter derivation for lattice-based e-voting: from a security target, the number of voters, the
number of authorities and the ballot format to the module rank, the modulus, the norm bounds, the
rejection parameters and the samplers.

The tool has two parts. A generic engine (`scripts/derive.py`) takes from a protocol module the M-SIS
bound, the feasibility conditions and the size as functions of the election and of the module rank,
searches every prime q = 17 (mod 32) below 2^50, and returns the smallest module rank and modulus that
reach separate hiding and binding targets under a chosen set of cost models. Two modules are
included: T-EVOLVE (threshold protocol, bounds from `param_report` of the T-EVOLVE code) and EVOLVE
(CCS 2017, constraints of its Section 5). This repository accompanies the paper "Deriving Lattice-Based
E-Voting Parameters for Different Election Sizes".

## Dependencies

| Dependency | Version | Used for |
|---|---|---|
| [t-evolve-cryptoalg](https://github.com/futures1mple/t-evolve-cryptoalg) | release v1.1, commit `1395c3f` (protocol code of release v1.0) | logged parameter chain, used as the baseline |
| [t-evolve-national-scale](https://github.com/futures1mple/t-evolve-national-scale) | commit `d0885a3` | `param_report`: bounds of the chain; identical to release v1.0 up to 2183 aggregation blocks (about 1.02 * 10^6 voters), where the tail exponent stays 187 ln 2 and two sampler levels suffice; adds the tail-exponent rule and a third sampler level for larger electorates |
| [lattice-estimator](https://github.com/malb/lattice-estimator) | commit `53da598` | hardness estimates of M-LWE and M-SIS (needs SageMath or passagemath) |

A C compiler, `make`, `git` and Python 3 are needed for the steps below.

## How to reproduce

The hardness estimates are stored in `results/estimates/` and `results/estimates_evolve/`, so every
step marked "no estimator" runs without Sage.

```sh
sh scripts/fetch_deps.sh               # clones the dependencies into deps/ and builds param_report
python3 scripts/check_baseline.py      # recomputes the logged configurations of release v1.1 (no estimator)
python3 scripts/sweep_electorate.py    # structural checks of the bound chain, about 10-15 minutes (no estimator)
```

Hardness tables (needs passagemath; see `SERVER_TASK.md` for running the full tables on a server):

```sh
python scripts/estimate_grid.py all --quick --jobs 2   # installation test, about 3 minutes
python scripts/estimate_grid.py all --jobs J           # integer-bit tables, about 4.5 hours on 7 cores
python scripts/evolve_chain.py estimate --jobs J       # M-SIS table of the EVOLVE module
python scripts/ceilings.py --jobs J                    # M-LWE at +-0.05 bit around each crossing
python scripts/derive.py estimate --jobs J             # M-SIS and M-LWE at the selected moduli
python scripts/naive_scaling.py --jobs J               # table of simpler strategies (2 M-LWE estimates)
python scripts/threshold_binding.py --jobs J           # binding and feasibility at the rank thresholds
```

Derivation, checks and figures (no estimator):

```sh
python3 scripts/check_monotone.py      # monotonicity of all stored estimates (hypothesis (H4))
python3 scripts/derive.py select       # results/derivation.{json,md}: selections, ceilings, thresholds
python3 scripts/figures.py             # results/figures/fig_scaling, fig_rank (.pdf, .png)
python3 scripts/paper_numbers.py      # results/paper_numbers.txt: table rows and derived numbers of the paper
```

## Scripts

| Script | Content |
|---|---|
| `derive.py` | the engine and the two protocol modules; targets: 128 bits under MATZOV, MATZOV with its quantum variant, core-SVP 0.292 beta, core-SVP 0.292 and 0.265 beta; 192 bits under MATZOV and MATZOV with quantum; hiding 192 / binding 128; ceilings with and without the dual attacks; thresholds by bisection |
| `estimate_grid.py` | M-LWE and M-SIS estimates at q_b, the largest prime q = 17 (mod 32) below 2^b |
| `ceilings.py` | hiding ceilings at 0.1-bit resolution |
| `check_monotone.py` | monotonicity of the estimates in q and in beta_SIS |
| `threshold_binding.py` | M-SIS estimate and the feasibility conditions at every rank threshold |
| `naive_scaling.py` | simpler strategies against the derived set at 6.48 million voters |
| `sweep_electorate.py` | bounds on a dense grid of 243 electorate sizes; monotonicity in N_V, d, n, L; independence of t and q |
| `check_baseline.py` | `param_report` against the logged configurations of release v1.1; expected output: `6 configurations reproduced, 0 differ` |
| `evolve_chain.py` | the EVOLVE constraints (used by the EVOLVE module of `derive.py`) and its M-SIS table; its `select` mode (moduli q_b only) is superseded by `derive.py` |
| `paper_numbers.py` | table rows and derived numbers of the paper from `derivation.json` and `param_report` |
| `figures.py` | Figure 1 of the paper (`fig_scaling`) and a plot of the selected rank (`fig_rank`) |
| `select_params.py` | superseded by `derive.py`: selection restricted to the moduli q_b (first manuscript version) |

## Layout

| Path | Content |
|---|---|
| `scripts/` | dependency setup, the engine, checks, hardness tables |
| `results/` | outputs of the scripts (see `results/README.md`) |
| `SERVER_TASK.md` | procedure for computing the hardness tables on a server |
| `deps/` | pinned dependencies, created by `fetch_deps.sh`, not tracked |

## License

MIT, see `LICENSE`.
