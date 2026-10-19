import numpy as np
import pylab as py


def load_spectrum(fname):
    x, f, c = np.loadtxt(fname, unpack=1, skiprows=1)
    lam = 3e10 / x * 1e8
    ff = f / lam**2
    s = np.argsort(lam)
    lam, ff = lam[s], ff[s]
    ff /= np.interp(1300.0, lam, ff)
    return lam, ff


def compare(pdf):
    """P-Cygni line test with resonance-line PRD (opacity_line_res_scatter=2).

    With epsilon=0 every line interaction is a scatter, so the
    equivalent width must match the coherent-scattering reference run
    (spectrum_coh.dat, res_scatter=1) to a few percent: PRD changes
    where the photon lands in the line, not the integrated flux. The
    wing/core flux ratio must also be physical -- wings populated but
    still depressed relative to the pseudo-continuum.
    """
    lam_p, ff_p = load_spectrum("spectrum_1.dat")
    lam_c, ff_c = load_spectrum("spectrum_coh.dat")

    xa, fa = np.loadtxt("analytic_line.dat", unpack=1)
    fa /= np.interp(1300.0, xa, fa)

    py.clf()
    py.plot(lam_p, ff_p)
    py.plot(lam_c, ff_c, "--")
    py.plot(xa, fa, linewidth=4, color="black")
    py.xlim(1050, 1350)
    py.ylim(0, 2)
    py.title("single line homologous, PRD vs coherent")
    py.xlabel("wavelength (angstroms)", size=15)
    py.ylabel("flux", size=15)
    py.legend(["SEDONA-PRD", "SEDONA-coherent", "analytic Sobolev"])
    if pdf != "":
        pdf.savefig()

    m = (lam_p > 1050) & (lam_p < 1350)
    mc = (lam_c > 1050) & (lam_c < 1350)
    ew_prd = np.trapezoid(1 - ff_p[m], lam_p[m])
    ew_coh = np.trapezoid(1 - ff_c[mc], lam_c[mc])

    failure = 0
    if abs(ew_prd - ew_coh) / abs(ew_coh) > 0.05:
        print(
            "FAIL: PRD EW %.1f A differs from coherent %.1f A by >5%%"
            % (ew_prd, ew_coh)
        )
        failure = 1

    wing = np.interp(1150.0, lam_p, ff_p)
    core = np.interp(1215.0, lam_p, ff_p)
    wing_core = wing / core
    if not 0.10 < wing_core < 0.45:
        print(
            "FAIL: wing/core flux ratio %.3f outside [0.10, 0.45]"
            % wing_core
        )
        failure = 1

    return failure


if __name__ == "__main__":
    status = compare("")
    print("SUCCESS" if status == 0 else "FAILURE")
