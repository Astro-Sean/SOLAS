import numpy as np
import pylab as py

M_P = 1.67262158e-24
K_B = 1.380658e-16
R0 = 1e3


def get_error(a, b, x=[], x_comp=[], use=[]):
    """Mean and max fractional error of a vs b (b interp'd to x)."""
    y_comp = np.interp(x, x_comp, b) if len(x) else b
    if len(use):
        a, y_comp = a[use], y_comp[use]
    err = np.abs(a - y_comp)
    return err.max() / y_comp.max(), err.mean() / y_comp.mean()


def exp_int_e3(tau):
    """2*E3(tau): slab transmission for cosine-law incident photons."""
    mu = np.linspace(1e-6, 1.0, 4000)
    return 2.0 * np.trapezoid(np.exp(-tau / mu) * mu, mu)


def compare(pdf, kappa=0.05):
    """Hydro state vs exact Riemann solution, then the spectrum.

    The transport stage emits a flat L_nu source at the inner boundary
    through the evolved shocked gas with a grey opacity. Photons are
    emitted cosine-law into the outward hemisphere, so each ray sees a
    slant column Sigma/mu and the transmitted fraction is 2*E3(kappa*Sigma),
    E3 being the third exponential integral.
    """
    failure = 0

    z, x, rho_e, p_e, v_e = np.loadtxt("exact_solution.dat", unpack=True)
    r, rho, v, T = np.loadtxt(
        "plt_00005.dat", unpack=True, usecols=[0, 1, 2, 3]
    )
    xme = r - R0
    p_me = rho / M_P * K_B * T
    t_e = p_e * M_P / rho_e / K_B
    use = xme < 5

    py.clf()
    py.suptitle("sod shock tube: hydro vs exact, spectrum through shock")
    for i, (ym, ye, lab) in enumerate(
        [
            (rho, rho_e, "density"),
            (v, v_e, "velocity"),
            (p_me, p_e, "pressure"),
            (T, t_e, "temperature"),
        ]
    ):
        py.subplot(2, 2, i + 1)
        py.xlim(0, 6)
        py.plot(x, ye, color="red", linewidth=3)
        py.plot(xme, ym, "o", markerfacecolor="none")
        py.ylabel(lab)
        max_err, mean_err = get_error(ym, ye, x=xme, x_comp=x, use=use)
        print("%s: mean_err=%.3f max_err=%.3f" % (lab, mean_err, max_err))
        if mean_err > 0.05:
            print("FAIL: %s mean error %.3f > 5%%" % (lab, mean_err))
            failure = i + 1
    if pdf != "":
        pdf.savefig()

    # emergent spectrum: measured transmission = attenuated/unattenuated
    nu, fl, cc = np.loadtxt("spectrum_1.dat", unpack=1, skiprows=1)
    nu_r, fl_r, cc_r = np.loadtxt("spectrum_ref.dat", unpack=1, skiprows=1)
    lam = 3e10 / nu * 1e8
    lam_r = 3e10 / nu_r * 1e8
    s = np.argsort(lam)
    lam, fl = lam[s], fl[s]
    fl_r = np.interp(lam, lam_r, fl_r)

    flam = fl * (nu[s] / lam) * 1e-8  # erg/s/angstrom
    flam_r = fl_r * (nu[s] / lam) * 1e-8

    r_emit = 990.3 * 1.01
    mg = r > r_emit
    sigma = np.trapezoid(rho[mg], r[mg])
    tau = kappa * sigma
    trans_pred = exp_int_e3(tau)
    inband = np.ones_like(lam, dtype=bool)
    trans_meas = np.median(fl[inband]) / np.median(fl_r[inband])
    print(
        "Sigma=%.3f tau=%.3f predicted trans=%.3f measured=%.3f"
        % (sigma, tau, trans_pred, trans_meas)
    )
    if abs(trans_meas - trans_pred) / trans_pred > 0.10:
        print(
            "FAIL: transmission %.3f vs analytic 2*E3(tau)=%.3f"
            % (trans_meas, trans_pred)
        )
        failure = 5

    py.clf()
    py.step(lam, flam, where="mid")
    py.step(
        lam,
        flam_r * trans_pred,
        color="black",
        linewidth=2,
        alpha=0.6,
        where="mid",
    )
    py.xlabel("wavelength (angstroms)", size=15)
    py.ylabel("L_lambda (erg/s/angstrom)", size=15)
    py.title(
        "spectrum through evolved shock, grey kappa=%.3f cm2/g" % kappa
    )
    py.legend(["SEDONA kappa>0", "2*E3(tau) analytic"])
    if pdf != "":
        pdf.savefig()

    return failure
