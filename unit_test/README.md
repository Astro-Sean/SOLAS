# unit_test

The single test suite for this code. Each unit test is a directory
containing `param.lua` and `run_test.py`. Every test keeps its own
output inside its folder: the run log and comparison figures go to
`<test>/plots/`.

## Layout

- `scripts/`                    runner and helper scripts
- `P-Cygni_line_prd/`           resonance-line PRD test: EW conservation
                                vs the coherent-scattering mode, wing/core
                                flux ratio, and offgrid=0
- `P-Cygni_line_prd/plots/`     that test's output figures

## Usage

From this directory:

    python scripts/run_unit_test.py P-Cygni_line_prd
    python scripts/run_unit_test.py P-Cygni_line_prd -n 4
    python scripts/run_unit_test.py --list

The script builds `mpirun -np <n> ../src/sedona6.ex param.lua`, runs it
inside the test directory, calls that test's `run_test.py`, and writes
figures to `P-Cygni_line_prd/plots/P-Cygni_line_prd.pdf`. The runner
always appends a page plotted from the raw output, so every run
produces a plot even if the test draws none itself.

`SEDONA_HOME` should point at the repository root, and `sedona6.ex`
should exist in `src/` (build with `cd src && ./install.sh general`).

## Adding a test

Make a new directory here with `param.lua` (the model input) and
`run_test.py` (exposing `run_test(pdf, runcommand)` that returns 0 on
success). Reference curves and model files belong in the test
directory; generated outputs (`plt_*`, `spectrum_*`, `plots/*.pdf`,
binaries) are not committed.
