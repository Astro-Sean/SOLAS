<img  align="right"  src="doc/sedona_logo.png" alt="drawing" width="200"/>

## Overview

SOLAS (Irish for "light") is a development of the SEDONA Monte Carlo
radiative transfer code by Kasen (https://github.com/dnkasen/pubsed).

Like SEDONA, SOLAS calculates the radiation signatures of supernovae
and other transient phenomena, on a personal computer or cluster, with
hydrodynamic solvers and an implicit Monte Carlo radiation transport
scheme.  See Kasen et al. (2006), Roth & Kasen (2015) for the method.

## New features relative to pubsed

Line physics:
- Resonance-line scattering channel: absorbed line photons can be
  re-emitted coherently instead of always thermalizing
  (`opacity_line_res_scatter`), with an optional partial-redistribution
  frequency shift (`opacity_line_prd_dnu`)
- Fluorescence channel: absorbed line photons can be re-emitted on
  other downward transitions of the same atom
- Branching ratios from van Regemorter collisional quenching versus
  total radiative decay rates
- Local Sobolev velocity gradient (|dv/dr|) for expansion-opacity
  line optical depths in non-homologous flows
- Configurable line wing width in Doppler units
  (`opacity_line_wing_dnu`)

Discrete diffusion (DDMC):
- Selectable stuck-packet handling (`ddmc_stuck_method`):
  teleport to a diffusion front (upstream behavior), hand back to IMC
  in place, or hold to census
- Optional multi-group DDMC with per-group Rosseland mean opacities
  (`ddmc_n_groups`)
- Rosseland-mean leakage opacity and DDMC diagnostics (entries,
  diffusion steps, handbacks, maximum steps)

Emission:
- Optional importance-weighted (flattened) radioactive emission
  (`particles_emit_flatten`)

Correctness and robustness:
- Energy-closure diagnostic per step (`debug_energy_closure`): emitted
  vs escaped vs absorbed-tally vs census energy
- Global MPI normalization of the escaped-spectrum rescale factor;
  previously each rank rescaled by its own noisy local fraction
- MPI reductions added for the J_nu and radiation-force tallies
- Per-thread Monte Carlo tallies remove OpenMP atomic contention
- Fail-safe temperature solver: failed root finds keep the previous
  gas temperature and emit a warning instead of silently zeroing it
- NLTE level-population solver falls back to LTE on singular rate
  matrices instead of aborting
- Inner-boundary velocity read from the model instead of an artificial
  zero (removes a spurious dv/dr ramp in truncated models)
- Comoving-frame frequency used in the DDMC interface escape check
- Out-of-bounds CDF sampling fixed; off-grid escaped packets counted
  and reported

## Get Started!

To get started, read the documentation at
https://dnkasen.github.io/sedona_documentation/
The documentation can be edited by modifying the .rst files the
doc/sphinx_doc folder

You can find (some of the) development guidelines at
https://github.com/dnkasen/sedona6/wiki/Developers-Guide.

To download additional atomic data, go to
https://www.dropbox.com/sh/pw9ttd90w503vzq/AAD3H78_5dGSFCRgBoupg0Iva?dl=0
