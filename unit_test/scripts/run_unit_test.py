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

    pdf.close()
    print("%s -> %s (status %d)" % (args.test, pdfname, status))
    return 0 if status == 0 else status


if __name__ == "__main__":
    sys.exit(main())
