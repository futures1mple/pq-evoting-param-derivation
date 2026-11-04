# Computing the hardness tables on a server

The integer-bit hardness tables (`results/estimates/lwe.jsonl`, 420 entries, and `sis.jsonl`, 1856
entries) take several hours. This is the procedure used for them, on a Linux server with 8 cores, in
the directory of this repository (`~/pq-evoting-param-derivation`). After each step, one line (date,
step, outcome) was appended to `results/estimates/PROGRESS.txt`.

Long runs should go inside a `tmux` session so that they survive a lost connection.

## Step 0. Machine record

    mkdir -p results/estimates
    { date -u; uname -a; nproc; free -g; df -h ~; uptime; python3 --version; } > results/estimates/machine.txt

Requirements: at least 8 GB free disk and 2 GB free memory per parallel process (see Step 4).

## Step 1. Dependencies and baseline

    sh scripts/fetch_deps.sh
    python3 scripts/check_baseline.py

`fetch_deps.sh` needs `git`, `make` and a C compiler (`build-essential`). `check_baseline.py` must end
with `6 configurations reproduced, 0 differ`.

## Step 2. Python environment with passagemath

Python 3.13 in a virtual environment at `deps/venv` (inside `deps/`, which is not tracked):

    uv venv --python 3.13 deps/venv
    . deps/venv/bin/activate
    uv pip install passagemath-standard==10.8.12
    python -c "import sage.all; print('sage ok')"

The installation takes about 4 GB. Without `uv`, `python3.13 -m venv deps/venv` and
`pip install passagemath-standard==10.8.12` are equivalent.

## Step 3. Quick test (about 3 minutes)

    . deps/venv/bin/activate
    python scripts/estimate_grid.py all --quick --jobs 2

Expected output (the same values as in the logs of T-EVOLVE):

    [1/1, ...] ('sis', 10000, 4, 1, 7, 42, 'MATZOV') min=144.55 (...)
    [1/1, ...] ('lwe', 7, 1, 45, 'MATZOV') min=128.86 (...)

## Step 4. Full run

Number of parallel processes `J`: the smaller of `nproc - 1` and (free memory in GB) / 2, at least 1.

    tmux new -d -s ja-run "cd ~/pq-evoting-param-derivation && . deps/venv/bin/activate && \
      python scripts/estimate_grid.py all --jobs J 2>&1 | tee -a results/estimates/run.log"

The script first computes the bounds for every configuration with `param_report` (10 to 20 minutes,
mostly for 6.48 * 10^6 and 10^7 voters), then the SIS part (1856 entries, mostly seconds each) and the
LWE part (420 entries of about 2 to 5 minutes each on one core). With `J = 7` the run took 4 h 25 min.
Progress: `tail -3 results/estimates/run.log`. After an interruption, the same command resumes;
finished entries are skipped.

## Step 5. Checks after the run

    wc -l results/estimates/sis.jsonl results/estimates/lwe.jsonl
    grep -c '"min": null' results/estimates/sis.jsonl results/estimates/lwe.jsonl

`lwe.jsonl` must have 420 lines and `sis.jsonl` 1856. `results/estimates/SUMMARY.txt` records start
and end time, `J`, line counts, the number of `null` entries, the slowest entries and any warnings.

## Step 6. Package

    tar czf ~/ja_estimates.tgz results/estimates results/estimates_quick
