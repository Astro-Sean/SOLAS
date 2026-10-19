#!/usr/bin/env python
"""Run a single unit test and save its comparison figure.

Usage:
    python run_unit_test.py <test_dir> [-n NRANKS]
    python run_unit_test.py --list

<test_dir> is a path relative to unit_test/, e.g. P-Cygni_line_prd.
A unit test is a directory containing param.lua and run_test.py.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("test", nargs="?", help="test dir relative to unit_test/")
    parser.add_argument("-n", dest="nproc", default="2", help="mpi ranks")
    parser.add_argument("--list", action="store_true", help="list tests")
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, "..", ".."))
    tests_dir = os.path.join(repo, "unit_test")
    plots_dir = os.path.join(tests_dir, "plots")

    if args.list:
        for root, dirs, files in os.walk(tests_dir):
            if "run_test.py" in files and "param.lua" in files:
                print(os.path.relpath(root, tests_dir))
        return 0

    if not args.test:
        parser.print_help()
        return 1

    test_dir = os.path.join(tests_dir, args.test)
    if not os.path.isdir(test_dir):
        print("no such test: %s" % args.test)
        return 1

    exe = os.path.join(repo, "src", "sedona6.ex")
    if not os.path.isfile(exe):
        print("missing executable: %s (build first)" % exe)
        return 1

    os.system("cp %s %s" % (exe, test_dir))
    os.chdir(test_dir)
    sys.path.insert(0, test_dir)

    runcommand = "mpirun -np %s ./sedona6.ex param.lua" % args.nproc
    pdfname = os.path.join(
        plots_dir, args.test.replace("/", "_") + ".pdf"
    )
    pdf = PdfPages(pdfname)

    try:
        import run_test

        status = run_test.run_test(pdf, runcommand)
    except Exception:
        import traceback

        traceback.print_exc()
        status = 200

    # every test gets an accompanying plot: if its own run_test/compare
    # produced nothing, fall back to plotting the raw output so the pdf
    # is never empty
    import glob

    import numpy as np
    import pylab as plt

    def spectrum_page():
        files = sorted(glob.glob("spectrum_*.dat"))
        if not files:
            return False
        plt.clf()
        for f in files:
            d = np.loadtxt(f, skiprows=1)
            if d.ndim == 2 and d.shape[0] > 0:
                plt.plot(d[:, 0], d[:, 1], label=f)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("frequency (Hz)")
        plt.ylabel("flux (erg/s/Hz)")
        plt.title(args.test + ": output spectrum")
        plt.legend(fontsize=6)
        pdf.savefig()
        return True

    def profile_page():
        files = sorted(glob.glob("plt_*.dat"))
        if not files:
            return False
        plt.clf()
        for f in files[:8]:
            d = np.loadtxt(f, skiprows=2)
            if d.ndim == 2 and d.shape[0] > 0 and d.shape[1] > 3:
                plt.plot(d[:, 0], d[:, 3], label=f)
        plt.xscale("log")
        plt.xlabel("radius (cm)")
        plt.ylabel("gas temperature (K)")
        plt.title(args.test + ": output profile")
        plt.legend(fontsize=6)
        pdf.savefig()
        return True

    own_pages = pdf.get_pagecount()
    if not spectrum_page():
        profile_page()

    pdf.close()

    if own_pages == 0:
        print("note: %s made no figure of its own; "
              "wrote a fallback output plot" % args.test)

    print("%s -> %s (status %d)" % (args.test, pdfname, status))
    return 0 if status == 0 else status


if __name__ == "__main__":
    sys.exit(main())
