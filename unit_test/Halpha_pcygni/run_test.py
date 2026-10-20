import os
import subprocess
import sys
import matplotlib
import compare

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SEDONA_HOME = os.environ["SEDONA_HOME"]
HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(SEDONA_HOME, "src", "sedona6.ex")


def main():
    os.chdir(HERE)
    subprocess.check_call(["python3", "make_mod.py"])
    for stale in ["spectrum_1.dat", "spectrum_1.h5"]:
        if os.path.exists(stale):
            os.remove(stale)

    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    subprocess.check_call(
        [
            "mpiexec",
            "--mca",
            "pml",
            "ob1",
            "--mca",
            "btl",
            "self,vader,tcp",
            "-np",
            str(nproc),
            EXE,
            "param.lua",
        ],
        stdout=open("run_output.log", "w"),
        stderr=subprocess.STDOUT,
    )
    if not os.path.exists("spectrum_1.dat"):
        print("FAIL: engine did not produce spectrum_1.dat")
        return 1

    os.makedirs("plots", exist_ok=True)
    d = np.loadtxt("spectrum_1.dat")
    nu, Lnu, clk = d[:, 0], d[:, 1], d[:, 2]
    lam = 2.9979e18 / nu
    Llam = Lnu * nu / lam
    o = np.argsort(lam)
    lam, Llam = lam[o], Llam[o]

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].step(lam, Llam, "k-", lw=1)
    ax[0].set_xlim(5300, 7800)
    ax[0].axvline(6563, color="r", ls="--", lw=1)
    ax[0].set_xlabel("lambda [A]")
    ax[0].set_ylabel("L_lambda [erg/s/A]")
    ax[1].step(lam, Llam, "k-", lw=1)
    ax[1].set_xlim(5900, 7300)
    ax[1].axvline(6563, color="r", ls="--", lw=1, label="H-alpha 6563 A")
    ax[1].set_xlabel("lambda [A]")
    ax[1].set_ylabel("L_lambda [erg/s/A]")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig("plots/Halpha_pcygni.pdf")

    ok = compare.main(d)
    print("STATUS = %d" % ok)
    return ok


if __name__ == "__main__":
    sys.exit(main())
