"""Homologous constant-density shell for the Sobolev P-Cygni test.

100 zones, v = r/t_exp from 2e7 to 2e9 cm/s, t_exp = 10 days,
rho = 1e-18 g/cm3, T = 6000 K. The photosphere sits at the core
boundary r_core = 1e9*t_exp so only the vphot..vmax shell scatters.
"""

nx = 100
vmax = 2e9
texp = 10 * 3600.0 * 24
temp = 6e3
rho = 1e-18

with open("constant_standard.mod", "w") as fout:
    fout.write("1D_sphere standard\n")
    fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
    fout.write("1.1\n")
    for i in range(nx):
        v = (i + 1.0) * vmax / nx
        r = v * texp
        fout.write("%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r, v, rho, temp))
