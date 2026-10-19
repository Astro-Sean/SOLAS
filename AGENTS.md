# SOLAS (sedona_source) working notes

## Build

MPCDF environment: intel/2024.0 + impi or openmpi, GSL and HDF5 from
the znver4 stack. Lua 5.2.1 is bundled in `src/external/lua-5.2.1`
(build it there first with `make`).

```
export SEDONA_HOME=<this dir>
export GSL_DIR=/mpcdf/soft/RHEL_9/packages/znver4/gsl/intel_2024.0-2024.0.2/2.6
export HDF5_DIR=/mpcdf/soft/RHEL_9/packages/znver4/hdf5/intel_2024.0-2024.0.2/1.14.1
export LUA_DIR=$SEDONA_HOME/src/external/lua-5.2.1
cd src && make        # produces src/sedona6.ex and src/build/gomc
```

At run time the GSL/HDF5/openmpi `lib` dirs must be in
`LD_LIBRARY_PATH`, and `SEDONA_HOME` must be set (params resolve
`defaults/sedona_defaults.lua` and `data/` through it). Use the
openmpi `mpiexec` matching `libmpi.so.40`.

## Tests

`unit_test/` is the only test suite; the old `tests/` tree was removed.
Each test lives in its own folder and produces its own output plus
`plots/<test>.pdf` inside that folder:

```
cd unit_test
python3 scripts/run_unit_test.py <test_name> -n <mpi_procs>
```

Tests assert `offgrid=0` and the energy-closure diagnostic where
applicable. Do not weaken a comparison to make a test pass.

## Plot conventions

- Spectra are plotted in CGS units: wavelength in Angstrom and flux as
  erg/s/cm^2/A (F_lambda) or erg/s/A (L_lambda) - state which on the
  axis label.
- `spectrum_*.dat` columns are nu [Hz], L_nu [erg/s/Hz], packet
  clicks; convert with `L_lambda = L_nu * c / lambda^2`.
- Every unit test must emit at least one figure page; the runner adds
  a raw-output fallback page when the test draws none.

## C/C++ conventions

Keep ASCII source. Physics fixes and test-suite changes go in separate
commits against `development`.
