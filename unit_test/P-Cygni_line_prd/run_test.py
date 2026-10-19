import os
import re

import compare


def run_test(pdf="", runcommand=""):

    # the harness appends ">> suite_log 2>&1"; strip that and capture a
    # local log so the offgrid assertion can read our own output
    log = "run_output.log"
    cmd = re.sub(r"\s*>>.*$", "", runcommand) if runcommand else \
        "./sedona6.ex param.lua"

    # coherent-scattering reference (res_scatter=1) for the EW
    # conservation check; a scattering line must give the same
    # integrated profile in both modes
    os.system(cmd.replace("param.lua", "param_coh.lua") + " > /dev/null 2>&1")
    os.system("cp spectrum_1.dat spectrum_coh.dat")

    os.system(cmd + " > " + log + " 2>&1")

    failure = compare.compare(pdf)

    txt = open(log).read()
    m = re.search(r"offgrid=(\S+)", txt)
    if m and float(m.group(1)) != 0.0:
        print("FAIL: offgrid escapes = %s" % m.group(1))
        failure = 1
    return failure
