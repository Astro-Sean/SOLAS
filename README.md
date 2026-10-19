## SOLAS

SOLAS is a fork of the SEDONA Monte Carlo radiative transfer code
(https://github.com/dnkasen/pubsed). It computes light curves and
spectra of supernovae and other transients using implicit Monte Carlo
transport, with optional hydrodynamics.

The code here is under development.

## Changes relative to pubsed

Line opacity:

- Resonance-line scattering: a line photon can be re-emitted at the
  same frequency instead of being thermalized
  (`opacity_line_res_scatter`, `opacity_line_prd_dnu`)
- Fluorescence: a line photon can come out on a different downward
  transition of the same atom
- The split between scattering, fluorescence, and absorption uses the
  collisional quenching rate versus the total radiative decay rate
- Expansion opacity uses the local |dv/dr| for the Sobolev length
  instead of assuming 1/t
- Line wing width is a parameter (`opacity_line_wing_dnu`)

Discrete diffusion:

- `ddmc_stuck_method` selects what to do when a packet takes too many
  diffusion steps: move it to a front, return it to IMC, or hold it
  to census
- `ddmc_n_groups` enables multi-group DDMC with per-group Rosseland
  means
- DDMC diagnostics are printed each step (entries, steps, handbacks)

Emission:

- `particles_emit_flatten` uses importance weighting for radioactive
  emission

Bug fixes and robustness:

- Energy-closure diagnostic per step (`debug_energy_closure`)
- Escaped-spectrum rescale now uses global MPI counters; each rank
  previously rescaled by its own local escape fraction
- Missing MPI reductions added for J_nu and radiation-force tallies
- Per-thread tallies replace OpenMP atomics
- A failed temperature solve keeps the previous temperature and warns
  instead of silently setting it to zero
- The NLTE solver falls back to LTE on a singular rate matrix
- The inner-boundary velocity is read from the model instead of being
  set to zero
- Off-grid escaped packets are counted and reported
- Fixed an out-of-bounds CDF sample that inflated escaped packet
  energies

## Get Started

Documentation: https://dnkasen.github.io/sedona_documentation/

Additional atomic data:
https://www.dropbox.com/sh/pw9ttd90w503vzq/AAD3H78_5dGSFCRgBoupg0Iva?dl=0

## References

- Kasen, Thomas & Nugent 2006, ApJ 651, 366
- Roth & Kasen 2015, ApJS 217, 9
