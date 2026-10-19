#ifndef _TRANSPORT_H
#define _TRANSPORT_H

#include <vector>
#include <list>

#include <string>

#include "particle.h"
#include "grid_general.h"
#include "cdf_array.h"
#include "locate_array.h"
#include "thread_RNG.h"
#include "spectrum_array.h"
#include "GasState.h"
#include "ParameterReader.h"
#include "VoigtProfile.h"
#include "sedona.h"
#include "h5utils.h"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


using std::vector;

/* class RelativityData */
/* { */
/*   double gamma; */
/*   double beta; */
/*   double V[3]; */
/*   double D[3]; */
/*   double v_dot_D; */
/* }; */





class transport
{

 private:

  // arrays of particles
  std::vector<particle> particles;
  std::vector<particle> particles_new; // For debugging checkpointing
  std::vector<particle> particles_escaped;
  std::vector<particle> particles_escaped_new;
  int max_total_particles;

  // gas class for opacities
  vector<GasState> gas_state_vec_;

  // pointer to parameter reader
  ParameterReader* params_;
  // atomic data
  std::string atomdata_file_;
  AtomicData* atomic_data_;

  std::string escaped_particle_filename_;
  int save_escaped_particles_;
  double maxn_escaped_particles_;

  // MPI stuff
  int MPI_nprocs;
  int MPI_myID;
  int my_zone_start_, my_zone_stop_;
  double *src_MPI_block, *dst_MPI_block;
  double *src_MPI_zones, *dst_MPI_zones;
#ifdef MPI_PARALLEL
  MPI_Datatype MPI_real;
#endif

  // simulation parameters
  double step_size_;
  int    steady_state;
  int    radiative_eq;
  int    first_step_;
  int    verbose;
  int    omit_scattering_;
  int    store_Jnu_;
  int    core_fix_luminosity_;
  double maximum_opacity_;
  int    last_iteration_;
  int    omit_composition_decay_;
  int    compton_scatter_photons_;
  double fleck_alpha_;
  int    solve_Tgas_with_updated_opacities_;
  int    set_Tgas_to_Trad_;
  int    fix_Tgas_during_transport_;

  int use_nlte_;
  double nlte_dens_min_;  // minimum density for NLTE solve; below this, fall back to LTE

  // DDMC diagnostics: counters accumulated over one transport step.
  // ddmc_entries_   : number of MC -> DDMC transitions (zone tau > threshold)
  // ddmc_leaps_     : total DDMC leakage/scatter/stay iterations executed
  // ddmc_to_imc_    : number of DDMC -> MC transitions (zone tau <= threshold)
  // ddmc_max_leaps_ : largest number of DDMC iterations for one particle
  // ddmc_stuck_     : particles needing >1e5 DDMC iterations (near-trapped)
  long long ddmc_entries_;
  long long ddmc_leaps_;
  long long ddmc_to_imc_;
  long long ddmc_max_leaps_;
  long long ddmc_stuck_;
  long long ddmc_max_leaps_allowed_ = 100000;
  // handling of packets that exceed ddmc_max_leaps_allowed_:
  //   0 = teleport to the diffusion front (legacy heuristic)
  //   1 = finish the step as an IMC packet in place (no_ddmc flag)
  //   2 = hold the packet until census (fate = stopped)
  int ddmc_stuck_method_ = 0;

  // print per-step energy-closure accounting (start/escaped/tallied/
  // dead/census) when nonzero; set by debug_energy_closure param
  int debug_energy_closure_ = 0;
  // cumulative packet-energy bookkeeping for the closure diagnostic:
  // emitted (start-of-step), escaped, absorbed-fate and tally totals
  double e_cum_start_ = 0, e_cum_esc_ = 0, e_cum_abs_ = 0,
         e_cum_dead_ = 0;


  // current time in simulation
  double t_now_;

  // inner boundary
  double L_core_, r_core_, T_core_, time_core_, core_frequency_;
  // emission spectrum from inner boundary
  cdf_array<double> core_emission_spectrum_;
  // emission distribution across zones
  cdf_array<double> zone_emission_cdf_;

  // emission point sources
  int use_pointsources_;
  vector<double> pointsource_x_, pointsource_y_, pointsource_z_;
  vector<double> pointsource_L_, pointsource_T_;
  cdf_array<double> pointsource_emission_cdf_;
  cdf_array<double> pointsource_emission_spectrum_;
  double pointsources_L_tot_;


  // For sampling Maxwell-Boltzmann distribution for Compton scatterirng
  cdf_array<double> mb_cdf_;
  double mb_dv;

  // minimum and maximum temperatures
  double temp_max_value_, temp_min_value_;

  // class to hold output spectrum
  spectrum_array optical_spectrum;
  spectrum_array optical_spectrum_new;
  spectrum_array gamma_spectrum;
  spectrum_array gamma_spectrum_new;

  // random number generator
  mutable thread_RNG rangen;

   // random number generator
  //gsl_rng *rangen;

  // Voigt profile class
  VoigtProfile voigt_profile_;

  // pointer to grid
  grid_general *grid;

  // the frequency grid for emissivity/opacity (Hz)
  locate_array nu_grid_;

  // boundary conditions
  int boundary_in_reflect_, boundary_out_reflect_;

  // array to weight the emissivity(size of nu_grid)
  vector<real>  emissivity_weight_;

  // the zone opacity/emissivity variables
  vector< cdf_array<OpacityType> >    emissivity_;
  // importance-sampled emission: fraction emit_flatten_ of new thermal
  // particles are drawn from a uniform-in-bin distribution instead of
  // the emissivity CDF, with p.e rescaled by emis_prob/sample_prob.
  // Optical line bins carry ~1e-7 of the zone emissivity weight, so
  // unbiased emission essentially never creates line photons even
  // though the gas emits them; flattening fixes the statistics, not
  // the physics (mean emitted energy spectrum is unchanged).
  double emit_flatten_ = 0;
  vector< cdf_array<OpacityType> >    emit_cdf_;
  vector< vector<double> >            emit_wt_;
  vector< vector<OpacityType> >       abs_opacity_;
  vector< vector<OpacityType> >       scat_opacity_;
  // resonance-line scattering opacity: the A/(A+C) albedo-weighted
  // part of bound-bound opacity. Photons interacting through this
  // channel scatter coherently in the comoving frame (keep nu,
  // isotropize direction) instead of thermalizing. Empty unless
  // opacity_line_res_scatter is set.
  vector< vector<OpacityType> >       line_res_opacity_;
  // fluorescence opacity: the part of line absorption where the upper
  // level decays through a different transition. The photon is re-emitted
  // in a branching line, so its frequency is resampled from the zone
  // bound-bound emissivity rather than kept or thermalized.
  vector< vector<OpacityType> >       line_fluor_opacity_;
  // per-zone CDF of bound-bound emissivity only; kernel for the
  // fluorescence frequency resample. Using the total emissivity here
  // would dump the photons into the UV-dominated continuum.
  vector< cdf_array<OpacityType> >    line_emis_cdf_;
  int use_line_res_scatter_ = 0;
  // rest frequency of the dominant resonance-contributing line in each
  // (zone, bin), set by the opacity pass. PRD redistribution samples
  // the outgoing detuning around this nu_0 (Hummer RII-style), not
  // around the incoming frequency.
  vector< vector<double> >            res_line_nu0_;
  // fractional frequency half-width of the PRD redistribution window
  // around nu_0. The window only limits which bins the sampling loop
  // visits; the actual support is the line profile itself, so values
  // larger than the local line width are harmless.
  double res_prd_dnu_ = 0.02;
  // recoherence probability for PRD: with this probability a resonantly
  // scattered photon is re-emitted exactly at the line rest frequency
  // instead of being drawn from the line profile (the coherent branch
  // of Hummer RII). 0 = complete redistribution over the profile.
  double res_prd_recoherence_ = 0.0;
  // counts get_opacity calls where the resonance+fluorescence fractions
  // had to be clamped because they exceeded the non-thermal share
  long long res_chan_clamps_ = 0;
  // trailing EMA of the steady-state escape rescale factor, so a single
  // noisy iteration does not jitter the output spectrum
  double rescale_fac_ema_ = 0;
  vector<OpacityType> planck_mean_opacity_;
  vector<OpacityType> rosseland_mean_opacity_;
  vector< vector<real> > J_nu_;
  vector<real> compton_opac;
  vector<real> photoion_opac;

  // the following are only resized and used if gas_state_.use_nlte_ is set to 1
  vector<real> bf_heating;
  vector<real> bf_cooling;
  vector<real> ff_heating;
  vector<real> ff_cooling;
  vector<real> coll_cooling;

  // discrete diffusion probabilities
  vector<real> ddmc_P_up_, ddmc_P_dn_;
  vector<real> ddmc_P_adv_;
  vector<real> ddmc_P_abs_;
  vector<real> ddmc_P_stay_;
  vector<real> ddmc_use_in_zone_;
  int use_ddmc_;
  double ddmc_tau_;
  locate_array randomwalk_x;
  vector<double> randomwalk_Pescape;

  // frequency-group-lumped DDMC: nu_to_group_ maps each transport
  // frequency bin to a DDMC group; rosseland_group_opacity_[zone][g]
  // holds the group Rosseland mean. Active only when
  // ddmc_n_groups_ > 1; otherwise the legacy grey-Rosseland path runs.
  int ddmc_n_groups_ = 1;
  vector<int> nu_to_group_;
  vector< vector<OpacityType> > rosseland_group_opacity_;

  // diffusion opacity governing DDMC membership and leakage for a
  // packet in frequency bin i_nu of the given zone. With group
  // lumping this is the packet's group Rosseland mean; with a single
  // group it falls back to the arguments' frequency-dependent sigma_i
  // (membership) or the zone Rosseland mean (leakage), preserving the
  // historical behavior.
  double ddmc_opacity(int zone, int i_nu, double sigma_i) const
  {
    if (ddmc_n_groups_ > 1)
      return rosseland_group_opacity_[zone][nu_to_group_[i_nu]];
    return sigma_i;
  }
  double ddmc_leak_opacity(int zone, int i_nu) const
  {
    if (ddmc_n_groups_ > 1)
      return rosseland_group_opacity_[zone][nu_to_group_[i_nu]];
    return rosseland_mean_opacity_[zone];
  }

  // the radiation quantities in the zone
  vector <real> e_rad;
  // line mean intensity
  vector< vector<real> > line_J_;
  double line_velocity_width_;

  // Per-thread radiation tally buffers indexed [thread][zone]. The
  // propagation loops write here instead of using omp atomic on the
  // shared zone fields, which serializes when many threads tally to
  // the same optically thick zone. Buffers are zeroed and merged into
  // the zone fields once per transport step. Allocated lazily in
  // step() since the zone count is known only after the grid exists.
  vector< vector<double> > tl_e_abs_;
  vector< vector<double> > tl_J_nu0_;
  vector< vector<double> > tl_fx_;
  vector< vector<double> > tl_fy_;
  vector< vector<double> > tl_fz_;
  vector< vector<double> > tl_fr_;

  int tl_tid() const {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
  }

  // Cache of solved LTE gas states, one entry per zone. For LTE
  // (use_nlte_ == 0) solve_state ignores J_nu and depends only on
  // (density, temperature, mass fractions), so when those inputs are
  // unchanged the previous solve can be restored exactly. This avoids
  // re-solving the Brent/level-population problem every iteration in
  // steady-iterate runs where zone T/rho are fixed. NLTE zones are
  // never cached because they depend on J_nu.
  vector<GasState> lte_state_cache_;
  vector<double>   lte_cache_T_;
  vector<double>   lte_cache_rho_;
  vector< vector<double> > lte_cache_X_;
  vector<char>     lte_cache_valid_;

  // setup functions
  void setup_core_emission();
  void setup_pointsource_emission();
  void read_pointsource_params(ParameterReader* par);

  // opacity functions
  int   get_opacity(particle&, double, double&, double&);
  int   get_opacity(particle&, double, double&, double&, double&, double&);
  void   set_opacity(double dt);
  double klein_nishina(double);
  double blackbody_nu(double T, double nu);
  void   reduce_opacities();

  // creation of particles functions
  void   emit_particles(double dt);
  void   emit_inner_source(double dt);
  void   emit_radioactive(double dt);
  void   emit_thermal(double dt);
  void   emit_heating_source(double dt);
  void   emit_from_pointsoures(double dt);
  void   create_isotropic_particle(int,PType,double,double,int importance=0);
  void   initialize_particles(int);
  // returns the energy weight the sampled frequency carries (1 unless
  // importance sampling flattened the emission CDF)
  double sample_photon_frequency(particle*,int importance=0);

  // special relativistic functions
  void   transform_comoving_to_lab(particle*);
  void   transform_lab_to_comoving(particle*);
  void   lorentz_transform(particle*, int);
  double dshift_comoving_to_lab(particle*);
  double dshift_lab_to_comoving(particle*);
  double do_dshift(particle*, int);

  // sampling Maxwell-Boltzmann distribution for Compton scatterirng
  void setup_MB_cdf(double, double, int);
  void sample_MB_vector(double, double*, double*);

  //propagation of particles functions
  ParticleFate propagate(particle &p, double tstop);
  ParticleFate propagate_monte_carlo(particle &p, double dt);
  ParticleFate discrete_diffuse_IMD(particle &p, double tstop);
  ParticleFate discrete_diffuse_DDMC(particle &p, double tstop);
  ParticleFate discrete_diffuse_RandomWalk(particle &p, double tstop);
  int move_across_DDMC_interface(particle &p, int, double, double);
  void setup_RandomWalk();
  void compute_diffusion_probabilities(double dt);
  void sample_dir_from_blackbody_surface(particle*);
  int clean_up_particle_vector();

  // scattering functions
  ParticleFate do_scatter(particle*, double, double = 0.0, double = 0.0);
  void compton_scatter(particle*);
  void compton_scatter_photon(particle*);
  void isotropic_scatter(particle*, int);

  // radiation quantities functions
  void wipe_radiation();
  void reduce_radiation(double);
  void reduce_Tgas();
  void reduce_n_elec();
  void reduce_Lthermal();

  // solve equilibrium temperature
  int solve_state_and_temperature(GasState*, int); // calls gas state solve from within interative solution for tempreature. For now, temperature solve is always based on radiative equilibrium
  void solve_eq_temperature();
  double rad_eq_function_LTE(GasState*, int,double,int, int &);
  double rad_eq_function_NLTE(GasState*, int,double,int, int &);
  double temp_brent_method(GasState*, int,int, int &);
  void   report_temp_solve_failure(int cell, double T_kept,
                                   const char* reason);

 public:

  int write_levels;


  void set_last_iteration_flag()
    {last_iteration_ = 1;}

  //--------------------------------
  // constructor and defaults
  //------------------------------

  transport()
  {
    time_core_ = 0;
  }

  // destructor
  ~transport();

  //----- functions ---------------

  // set things up
  void init(ParameterReader*, grid_general*);
  void setup_MPI();

  // run a transport step
  void step(double dt);

  // return
  int n_particles() { return particles.size(); }
  int n_particles_new() { return particles_new.size(); }
  int n_particles_escaped() { return particles_escaped.size();}
  int n_particles_escaped_new() { return particles_escaped_new.size();}

  // finalize and output spectra
  void output_spectrum();
  void output_spectrum(int);

  // print out functions
  void write_levels_to_plotfile(int);
  void write_radiation_file(int);
  void wipe_spectra();
  void clearEscapedParticles();

  void writeCheckpointParticlesAll(std::string fname);
  void writeCheckpointParticles(std::vector<particle>& particle_list,
      std::string fname, std::string groupname);
  void writeParticleProp(std::string fname, std::string fieldname,
      std::string groupname, std::vector<particle>& particle_list,
      int total_particles, int offset);
  void writeCheckpointSpectra(std::string fname);
  void writeCheckpointRNG(std::string fname);

  void readCheckpointParticles(std::vector<particle>& particle_list, 
      std::string fname, std::string groupname, bool test=false,
      bool all_one_rank=false);
  void readParticleProp(std::string fname, std::string fieldname,
      std::string groupname, std::vector<particle>& particle_list,
      int total_particles, int offset);
  void readCheckpointSpectra(std::string fname, bool test=false);
  void readCheckpointRNG(std::string fname, bool test=false);

  void testCheckpointParticles(std::string fname);
  void testCheckpointSpectrum(std::string fname);
};


#endif
