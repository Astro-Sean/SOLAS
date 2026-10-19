#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import time
import optparse
from importlib import import_module
import timeit


###############################################
# runs a series of tests of sedona
# by default, the tests to be run are read
# from the file "suite_test_list"
#
# output figures showing test results are put
# into in the pdf file:
#  test_results_%date.pdf
# (where %date is the current date)
#
# Usage:
#  python run_test_suite.py [options]
# Options:
#
#   -n 2
#   (will run on 2 mpi ranks, more generally -n nranks
#   will run on nranks ranks where nranks is an integer)
#
#  --testlist "spherical_lightbulb/1D","lucy_supernova/1D"
#   (will run the tests given after the --testlist flag)
#
#  --testfile "my_tests.txt"
#   (will read the file "my_tests.txt" or whatever filename you wante
#   to get the list of files to run)
###########################################

parser = optparse.OptionParser()
parser.add_option("-n",dest="nproc")
parser.add_option("--testlist","-t",dest="testlist",type='string')
parser.add_option("--testfile","-f",dest="testfile",type="string")
parser.add_option("--outfile","-o",dest="outfile",type="string")
parser.add_option("-v",action="store_true",dest="verbose")

(opts, args) = parser.parse_args()

nproc = 1
if (opts.nproc): nproc = opts.nproc

## executable directory and file
exec_dir   = "../src/"
executable = "sedona6.ex"
paramname = "param.lua"

homedir = os.getcwd()
date = time.strftime("%m-%d-%y")

# setup output files
if (opts.outfile):
    outfile = homedir + '/' + opts.outfile + '.txt'
    pdffile = homedir + '/' + opts.outfile + '.pdf'
else:
    outfile = homedir + '/test_results_' + date + '.txt'
    pdffile = homedir + '/test_results_' + date + '.pdf'
pdf = PdfPages(pdffile)

runcommand = "mpirun -np " + str(nproc) + " ./" + executable + " " + paramname
if (not opts.verbose):
    runcommand = runcommand + " >> " + outfile + " 2>&1"


########################################

# get list of tests eihter from input flag
if (opts.testlist):
    testlist = opts.testlist.split(',')
# or else from the default list file
else:
    testfile_name = "suite_test_list"
    if (opts.testfile): testfile_name = opts.testfile
    fin = open(testfile_name,"r")
    testlist = []
    for line in fin:
        line = line.rstrip('\n')
        line = line.rstrip(' ')
        if (os.path.isdir(line)):
            testlist.append(line)

line = "echo \"Testing SEDONA code on " + date + "\"  > " + outfile
print("Testing SEDONA code on " + date + "\n")
os.system(line)

# get number of threads
try:
    nthreads = os.environ['OMP_NUM_THREADS']
except:
    nthreads = 1

print("Will run " + str(len(testlist)) + " tests on " +str(nproc)  + " mpi ranks,"),
print("with " + str(nthreads) + " threads per rank\n");
print("------------------------------------------")
i = 1
for this_test in testlist:
    print(str(i) +') ' + this_test)
    i = i + 1
print("------------------------------------------")
print("\n")


total_status = 0
cnt = 1
for this_test in testlist:

    starttime = timeit.default_timer()

    hdr = "\n\n------------------------------------\n"
    hdr = hdr + "-- RUNNING: " + this_test + "\n"
    hdr = hdr + "------------------------------------\n"
    cmd = "echo \"" + hdr + "\"" +  " >> " + outfile
    os.system(cmd)

    print("------------------------------------------")
    print("- test " + str(cnt) + ") " + this_test)

    # copy the executable in by fixed name and fail fast: a missing or
    # failed cp previously went unnoticed and the test ran on a stale
    # or absent binary
    exec_src = os.path.join(exec_dir, executable)
    exec_dst = os.path.join(this_test, executable)
    if not os.path.isdir(this_test):
      print("SKIPPED: no such test directory", this_test)
      total_status += 404
      continue
    cp_status = os.system("cp " + exec_src + " " + exec_dst)
    if cp_status != 0 or not os.path.isfile(exec_dst):
      print("CRASHED: could not copy executable", exec_src)
      os.chdir(homedir)
      total_status += 200
      cnt += 1
      continue
    os.chdir(this_test)

    # immport and run script; the import itself can throw (e.g. py2
    # syntax), so it must live inside the try
    sys.path.append(os.getcwd())
    try:
      import run_test
      status = run_test.run_test(pdf,runcommand)
    except Exception:
      import traceback
      traceback.print_exc()
      status = 200
    sys.path.remove(os.getcwd())
    # evict test-local helper modules too, or the next test's
    # "import compare" resolves to a stale module from a previous
    # test directory (compare.py names collide across tests)
    for m in ['run_test', 'compare']:
      if m in sys.modules:
        del sys.modules[m]

    stoptime = timeit.default_timer()
    print("time = {:.1f} seconds".format(stoptime - starttime))

    if (status == 0):
        print("PASSED")
    elif (status == 200):
        print("CRASHED")
    else:
        print("FAILED: error code = ", status)
    print("------------------------------------------\n")
    total_status += status

    # return home
    os.chdir(homedir)
    cnt += 1


pdf.close()
