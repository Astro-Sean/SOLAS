import os
import re

import numpy as np

import compare


def run_test(pdf="", runcommand=""):

    base = (
        re.sub(r"\s*>>.*$", "", runcommand)
        if runcommand
        else "./sedona6.ex param.lua"
    )

    failure = 0

    # stage 1: pure hydro Sod tube to t=0.5; the t=0.4 snapshot
    # (plt_00005.dat) is compared against the exact Riemann solution
    os.system(base.replace("param.lua", "param_hydro.lua") + " > hydro.log 2>&1")
    if not os.path.exists("plt_00005.dat"):
        print("FAIL: hydro stage produced no plot files")
        return 1

    # freeze the t=0.4 state into a transport model
    r, rho, v, T = np.loadtxt("plt_00005.dat", unpack=True, usecols=[0, 1, 2, 3])
    with open("evolved.mod", "w") as f:
        f.write("1D_sphere standard\n%d\t%e\t%e 1 \n1.1\n" % (len(r), r.min(), 0.0))
        for i in range(len(r)):
            f.write(
                "%16.10e %15.8e %16.8e %16.8e 1.0\n" % (r[i], v[i], rho[i], T[i])
            )

    # stage 2: flat source through the shocked gas, grey kappa; the
    # unattenuated reference (kappa=0) gives the transmission numerator
    for stale in ("spectrum_1.dat", "spectrum_1.h5"):
        if os.path.exists(stale):
            os.remove(stale)
    os.system(base.replace("param.lua", "param_spec.lua") + " > spec.log 2>&1")
    os.system("cp spectrum_1.dat spectrum_opac.dat")
    os.system(base.replace("param.lua", "param_ref.lua") + " > spec_ref.log 2>&1")
    os.system("cp spectrum_1.dat spectrum_ref.dat")
    os.system("cp spectrum_opac.dat spectrum_1.dat")

    txt = open("spec.log").read()
    if "CALCULATION took" not in txt:
        print("FAIL: transport run did not complete")
        return 1

    failure = compare.compare(pdf)

    m = re.search(r"offgrid=(\S+)", txt)
    if m and float(m.group(1)) != 0.0:
        print("FAIL: offgrid escapes = %s" % m.group(1))
        failure = 1
    return failure
