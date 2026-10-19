import numpy as np
import pylab as py

import make_analytic


def load_spectrum(fname):
    x, f, c = np.loadtxt(fname, unpack=1, skiprows=1)
    lam = 3e10 / x * 1e8
    s = np.argsort(lam)
    return lam[s], x[s], f[s]


def compare(pdf):
    """Compare the emergent P-Cygni profile to the Sobolev analytic.

    The model is a single resonance line (2-level atom, 1215.7 A) in a
    homologous shell at effective tau0 ~ 0.3, i.e. the single-flight
    regime where the escape-probability analytic is valid. Shape
    metrics are evaluated on the F_nu profile (the convention the
    analytic integral shares); the figure itself is plotted in CGS:
    wavelength in angstroms and L_lambda in erg/s/angstrom.
    """
    lam, nu, fnu = load_spectrum("spectrum_1.dat")
    cont = np.median(fnu[(lam > 1050) & (lam < 1110)])
    fn = fnu / cont

    xa = np.arange(1050.0, 1400.0, 1.0)
    fa = np.array([make_analytic.sobolev_flux(l) for l in xa])

    # L_lambda [erg/s/angstrom] = L_nu * nu / lambda for plotting;
    # the analytic curve is converted with the same factor
    flam = fnu * (nu / lam) * 1e-8
    xa_nu = 3e10 / xa * 1e8
    fa_lam = fa * cont * (xa_nu / xa) * 1e-8

    py.clf()
    py.plot(lam, flam)
    py.plot(xa, fa_lam, linewidth=3, color="black", alpha=0.6)
    py.xlim(1080, 1350)
    py.title("single line homologous, Sobolev analytic comparison")
    py.xlabel("wavelength (angstroms)", size=15)
    py.ylabel("L_lambda (erg/s/angstrom)", size=15)
    py.legend(["SEDONA", "analytic Sobolev"])
    if pdf != "":
        pdf.savefig()

    failure = 0
    m = (lam > 1134) & (lam < 1298)
    fa_i = np.interp(lam[m], xa, fa)
    rms = np.sqrt(np.mean((fn[m] - fa_i) ** 2))
    print("profile rms vs analytic = %.4f" % rms)
    if rms > 0.05:
        print("FAIL: profile rms %.4f exceeds 0.05" % rms)
        failure = 1

    trough_lam = lam[m][fn[m].argmin()]
    trough = fn[m].min()
    print("trough %.3f at %.0f A" % (trough, trough_lam))
    if not (1130 < trough_lam < 1170):
        print("FAIL: absorption trough at %.0f A, expected ~1146" % trough_lam)
        failure = 1
    if not (0.4 < trough < 0.95):
        print("FAIL: trough depth %.3f outside [0.4, 0.95]" % trough)
        failure = 1

    mb = m & (lam > 1205) & (lam < 1290)
    bump_lam = lam[mb][fn[mb].argmax()]
    bump = fn[mb].max()
    print("bump %.3f at %.0f A" % (bump, bump_lam))
    if not (1205 < bump_lam < 1265):
        print("FAIL: emission bump at %.0f A, expected ~1225" % bump_lam)
        failure = 1
    if not (1.0 < bump < 1.5):
        print("FAIL: bump height %.3f outside [1.0, 1.5]" % bump)
        failure = 1

    return failure


if __name__ == "__main__":
    status = compare("")
    print("SUCCESS" if status == 0 else "FAILURE")
