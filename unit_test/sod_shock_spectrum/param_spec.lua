sedona_home = os.getenv("SEDONA_HOME")
defaults_file = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file = sedona_home.."/data/2level_atomdata.hdf5"
grid_type  = "grid_1D_sphere"
model_file = "evolved.mod"

days = 3600.0*24
core_n_emit      = 1e6
core_radius      = 990.3
core_luminosity  = 1e42
core_temperature = 0
core_spectrum_file = "core_spec.dat"

transport_nu_grid   = {2e15,6e15,0.001,1}
spectrum_nu_grid   = {2e15,6e15,4e12}

opacity_epsilon             = 1.0
opacity_grey_opacity        = 0.05
opacity_bound_bound         = 0
opacity_line_expansion      = 0

particles_n_emit_thermal    = 0
particles_n_emit_radioactive = 0

transport_steady_iterate        = 1
transport_radiative_equilibrium = 0
