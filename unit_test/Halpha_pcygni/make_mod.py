"""Homologous shell for the H-alpha P-Cygni demo.

100 zones, v = r/t_exp from 2e7 to 2e9 cm/s, t_exp = 10 days,
rho = 5e-13 g/cm3, T = 6000 K. Photosphere at the core boundary
r_core = 1e9*t_exp. Density is chosen so hydrogen is partially
neutral at this T, giving a saturated H-alpha line.
"""

nx = 100
vmax = 2e9
texp = 10 * 3600.0 * 24
temp = 6e3
rho = 5e-13

with open("constant_standard.mod", "w") as fout:
    fout.write("1D_sphere standard\n")
    fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
    fout.write("1.1\n")
    for i in range(nx):
        v = (i + 1.0) * vmax / nx
        r = v * texp
        fout.write("%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r, v, rho, temp))
