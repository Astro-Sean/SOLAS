# unit_test

The single test suite for this code. Each unit test is a directory
containing `param.lua` and `run_test.py`; running a test produces its
comparison figures in `plots/`.

## Layout

- `scripts/`            runner and helper scripts
- `plots/`              output figures land here (not tracked)
- `P-Cygni_line_prd/`   resonance-line PRD test: EW conservation vs the
                        coherent-scattering mode, wing/core flux ratio,
                        and offgrid=0

## Usage

From this directory:

    python scripts/run_unit_test.py P-Cygni_line_prd
    python scripts/run_unit_test.py P-Cygni_line_prd -n 4
    python scripts/run_unit_test.py --list

The script builds `mpirun -np <n> ../src/sedona6.ex param.lua`, runs it
inside the test directory, calls that test's `run_test.py`, and writes
the comparison figure to `plots/<test>.pdf`.

`SEDONA_HOME` should point at the repository root, and `sedona6.ex`
should exist in `src/` (build with `cd src && ./install.sh general`).

## Adding a test

Make a new directory here with `param.lua` (the model input) and
`run_test.py` (exposing `run_test(pdf, runcommand)` that returns 0 on
success). Reference curves and model files belong in the test directory;
generated outputs (`plt_*`, `spectrum_*`, binaries) are not committed.
