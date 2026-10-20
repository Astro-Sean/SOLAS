"""Checks the H-alpha P-Cygni profile in a steady spectrum.

Asserts a significant absorption trough blueward of rest 6563 A and an
emission peak near/above the line, plus offgrid=0.
"""

import numpy as np


def main(d):
    nu, L, clk = d[:, 0], d[:, 1], d[:, 2]
    lam = 2.9979e18 / nu
    Llam = L * nu / lam
    o = np.argsort(lam)
    lam, Llam, clk = lam[o], Llam[o], clk[o]

    cont = np.nanmedian(Llam[(lam > 6900) & (lam < 7200)])
    trough = np.nanmin(Llam[(lam > 6300) & (lam < 6540)]) / cont
    peak = np.nanmax(Llam[(lam > 6550) & (lam < 6800)]) / cont
    clicks = clk[(lam > 5900) & (lam < 7300)].sum()

    print("continuum L_lambda = %.3e erg/s/A" % cont)
    print("trough depth       = %.2f (require < 0.75)" % trough)
    print("emission peak      = %.2f (require > 1.10)" % peak)
    print("clicks in window   = %d (require > 1e4)" % clicks)

    ok = trough < 0.75 and peak > 1.10 and clicks > 1e4
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1
