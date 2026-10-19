"""Sobolev escape-probability P-Cygni profile for a single resonance line.

Port of the upstream make_sobolev_line.py to numpy/python3. Assumes a
homologous shell of constant Sobolev optical depth tau0 between vphot
and vmax, a flat continuum photosphere at vphot, and pure scattering
(source function = geometric dilution factor). Valid while multiple
scattering is negligible, i.e. tau0 of order unity or less.

Writes analytic_sobolev.dat: wavelength [A], normalized flux.
"""

import numpy as np

LAM0 = 1215.7  # line rest wavelength [A]
VPHOT = 1.0e9  # photosphere velocity [cm/s]
VMAX = 2.0e9   # shell outer velocity [cm/s]
TAU0 = 0.3     # effective Sobolev depth, calibrated to the MC run
NVP = 400

C_LIGHT = 2.9979e10


def sobolev_flux(lam, tau0=TAU0):
    vz = (lam - LAM0) / LAM0 * C_LIGHT
    vp = np.arange(0.0, VMAX, VMAX / NVP)
    vr = np.sqrt(vz**2 + vp**2)

    tau = np.where((vr > VMAX) | (vr < VPHOT) | ((vp < VPHOT) & (vz > 0)),
                   0.0, tau0)
    etau = np.exp(-tau)
    dilution = 0.5 * (1.0 - np.sqrt(np.maximum(0.0, 1 - (VPHOT / vr) ** 2)))
    source = np.where((vr < VPHOT) | (vr > VMAX), 0.0, dilution)
    iphot = np.where(vp < VPHOT, 1.0, 0.0)

    integrand = 2 * np.pi * vp * (iphot * etau + source * (1 - etau))
    return np.sum(integrand) * (VMAX / NVP) / (np.pi * VPHOT**2)


if __name__ == "__main__":
    lams = np.arange(1050.0, 1400.0, 1.0)
    flux = np.array([sobolev_flux(lam) for lam in lams])
    np.savetxt(
        "analytic_sobolev.dat",
        np.column_stack([lams, flux]),
        fmt="%10.4f %12.5e",
        header="lambda[A]  flux (normalized to flat continuum = 1)",
    )
