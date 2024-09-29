<p align="center">
  <img src="./logo.png" alt="SOLAS Logo" width="400">
</p>

<div align="center">
  <h4>The Concept</h4>
</div>

**SOLAS** (meaning *Light* in Irish) is an open-source ray tracing code developed to analyze emission line profiles in transient astronomical events.

We are entering a new era of supernova science, where the progenitor may be spectroscopically observed before it undergoes core-collapse (see [SN 2023fyq](https://arxiv.org/abs/2401.15148)). Although an exciting concept, the models needed to understand these observations are sparse. Motivated by this, *SOLAS* is conceptualized to provide a solution.

Currently in development, *SOLAS* will soon release a minimal viable product, building on the legacy of [*escatter.py*](https://github.com/Astro-Sean/escatter).

<div align="center">
  <h4>Challenges in Transient Astronomy</h4>
</div>

Most existing tools for studying transient events assume spherical symmetry, which fails to represent the highly asymmetric nature of these explosions. While several codes can model these phenomena, **SOLAS** aspires to be the first open-source, user-friendly tool that effectively incorporates the complexities of asymmetry in transient events.

<div align="center">
  <h4>Features Under Development</h4>
</div>

- Systems will comprise emission from H&alpha; 6563&Aring; and/or He I 5875&lambda; photons.
- Trace a photon from its emission at the photosphere through an intervening medium, capturing scattering, emission, and absorption events, as well as occultation effects.
- Generate detailed emission and absorption line profiles to analyze the intervening material and, crucially, its geometry.
- Position an "observer" at various orientations to assess line-of-sight effects.
- Utilize your own models — if you have data on the location of your photosphere, along with the density, distribution, and temperature of the intervening material — to explore potential emergent profiles.

<div align="center">
  <h4>Get Involved</h4>
</div>

If you’re interested in contributing to the development of this open-source tool, please [contact me](mailto:sean.brennan@astro.su.se).
