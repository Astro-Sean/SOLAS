"""Checks the H-alpha emission+absorption profile in a steady spectrum.

The setup combines a core continuum with thermal envelope emission, so
the profile should show a strong emission line near 6563 A with a
blueshifted absorption notch on its shoulder, plus offgrid=0.
"""

import numpy as np


def main(d):
    nu, L, clk = d[:, 0], d[:, 1], d[:, 2]
    lam = 2.9979e18 / nu
    Llam = L * nu / lam
    o = np.argsort(lam)
    lam, Llam, clk = lam[o], Llam[o], clk[o]

    cont = np.nanmedian(Llam[(lam > 6900) & (lam < 7200)])
    trough = np.nanmin(Llam[(lam > 6280) & (lam < 6460)]) / cont
    peak = np.nanmax(Llam[(lam > 6550) & (lam < 6700)]) / cont
    clicks = clk[(lam > 5900) & (lam < 7300)].sum()

    print("continuum L_lambda = %.3e erg/s/A" % cont)
    print("absorption notch   = %.2f (require < 0.55)" % trough)
    print("emission peak      = %.2f (require > 1.15)" % peak)
    print("clicks in window   = %d (require > 1e4)" % clicks)

    ok = trough < 0.55 and peak > 1.15 and clicks > 1e4
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1
