# unit_test

Lightweight harness for running a single test from `tests/` and saving
its comparison figures here, without running the whole suite.

## Layout

- `scripts/`  runner and helper scripts
- `plots/`    output figures land here (not tracked)

## Usage

From this directory:

    python scripts/run_unit_test.py P-Cygni_line
    python scripts/run_unit_test.py adiabatic_expansion/ddmc -n 4
    python scripts/run_unit_test.py --list

The script builds `mpirun -np <n> ../src/sedona6.ex param.lua`, runs it
inside the chosen test directory, calls that test's `run_test.py`, and
writes the comparison figure to `plots/<test>.pdf`.

`SEDONA_HOME` should point at the repository root, and `sedona6.ex`
should exist in `src/` (build with `cd src && ./install.sh general`).
