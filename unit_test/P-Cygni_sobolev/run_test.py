import os
import re

import compare


def run_test(pdf="", runcommand=""):

    log = "run_output.log"
    cmd = (
        re.sub(r"\s*>>.*$", "", runcommand)
        if runcommand
        else "./sedona6.ex param.lua"
    )

    # remove stale outputs so compare cannot read a previous run's
    # spectrum if this engine run fails or is interrupted
    for stale in ("spectrum_1.dat", "spectrum_1.h5"):
        if os.path.exists(stale):
            os.remove(stale)

    os.system(cmd + " > " + log + " 2>&1")

    txt = open(log).read()
    failure = 0
    if "CALCULATION took" not in txt or not os.path.exists("spectrum_1.dat"):
        print("FAIL: engine run did not complete")
        return 1

    failure = compare.compare(pdf)

    m = re.search(r"offgrid=(\S+)", txt)
    if m and float(m.group(1)) != 0.0:
        print("FAIL: offgrid escapes = %s" % m.group(1))
        failure = 1

    m = re.search(r"residual=(-?\S+) \((-?\S+) frac\)", txt)
    if m and abs(float(m.group(2))) > 1e-3:
        print("FAIL: energy closure residual frac = %s" % m.group(2))
        failure = 1
    return failure
