"""Sod shock tube initial model (Riemann problem at r = 2 cm).

500 zones spanning r = 1e3..1e3+10 cm; the large inner offset keeps
spherical curvature negligible so the planar exact solution applies.
Left state: rho=10, P=100. Right state: rho=1, P=1.
"""

nx = 500
rho = 1
r0 = 1e3
rsh = 10
dr = rsh / (1.0 * nx)
v = 0

k = 1.380658e-16
m_p = 1.67262158e-24

with open("sod_model.mod", "w") as fout:
    fout.write("1D_sphere standard\n")
    fout.write(str(nx) + "\t" + str(r0) + "\t" + str(0.0) + " 1 \n")
    fout.write("1.1\n")
    for i in range(nx):
        r = (i + 1) * dr
        if r < 2:
            rho = 10.0
            p = 100.0
        else:
            rho = 1.0
            p = 1.0
        T = p / k / (rho / m_p)
        fout.write("%16.10e %15.8e %16.8e %16.8e 1.0\n" % (r0 + r, v, rho, T))
