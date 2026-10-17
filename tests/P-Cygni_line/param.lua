sedona_home = os.getenv("SEDONA_HOME")
defaults_file = sedona_home.."/defaults/sedona_defaults.lua"
-- sedona_defaults.lua points at data/cmfgen_atomdata.hdf5, which is not
-- shipped; override with the bundled atom_data.hdf5
data_atomic_file = sedona_home.."/data/atom_data.hdf5"
grid_type  = "grid_1D_sphere"           -- grid geometry; must match input model
model_file = "constant_standard.mod"    -- input model file

-- helper variables
days = 3600.0*24

-- inner source
core_n_emit      = 2e6
core_radius      = 1e9*10*days
core_luminosity  = 1e42
core_temperature = 0

-- frequency grid
nu2 = 3e10/(1050*1e-8)
nu1 = 3e10/(1400*1e-8)
dnu = (nu2-nu1)/500.0
transport_nu_grid   = {nu1,nu2,0.0002,1}
spectrum_nu_grid   = {nu1, nu2*1.1, dnu}

opacity_epsilon             = 0.0
opacity_line_expansion      = 0
opacity_bound_bound         = 1
line_velocity_width         = 1e7

transport_steady_iterate        = 1
transport_radiative_equilibrium = 0



