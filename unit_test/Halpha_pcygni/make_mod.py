"""Homologous power-law shell for the H-alpha P-Cygni demo.

200 zones, v = r/t_exp from 1e7 to 2e9 cm/s, t_exp = 10 days.
rho(v) = rho0 * (v/v_ref)^-n, n = 10, as in RSG envelope models.
T = 6000 K uniform. The photosphere sits at the core boundary
r_core = 1e9*t_exp, so only the vphot..vmax shell forms the line.
"""

nx = 200
vmax = 1.3e9
vref = 1e9
texp = 10 * 3600.0 * 24
temp = 6e3
rho0 = 1.5e-11
npow = 10.0

with open("constant_standard.mod", "w") as fout:
    fout.write("1D_sphere standard\n")
    fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
    fout.write("1.1\n")
    for i in range(nx):
        v = (i + 1.0) * vmax / nx
        r = v * texp
        rho = rho0 * max(v / vref, 1.0) ** (-npow)
        fout.write("%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r, v, rho, temp))
