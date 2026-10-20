sedona_home = os.getenv("SEDONA_HOME")
defaults_file = sedona_home.."/defaults/sedona_defaults.lua"
-- 3-level hydrogen atom: Ly-alpha, Ly-beta and H-alpha (n2->n3)
data_atomic_file = sedona_home.."/data/H_3lev_plusC_atomdata.hdf5"
grid_type  = "grid_1D_sphere"
model_file = "constant_standard.mod"

days = 3600.0*24

core_n_emit      = 1e6
core_radius      = 1e9*10*days
core_luminosity  = 3e44
core_temperature = 0
core_spectrum_file = "core_spec.dat"

-- optical window around H-alpha (6563 A)
nu1 = 3e10/(7800*1e-8)
nu2 = 3e10/(5300*1e-8)
dnu = (nu2-nu1)/1000.0
transport_nu_grid = {nu1,nu2,0.0002,1}
spectrum_nu_grid  = {nu1, nu2*1.05, dnu}

opacity_epsilon          = 0.0
opacity_line_expansion   = 0
opacity_bound_bound      = 1
opacity_line_res_scatter = 1
line_velocity_width      = 1e7

transport_steady_iterate        = 1
transport_radiative_equilibrium = 0
particles_n_emit_thermal = 5e5
particles_emit_flatten = 0.3
