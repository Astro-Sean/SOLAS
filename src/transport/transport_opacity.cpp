#include <math.h>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <map>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "transport.h"
#include "physical_constants.h"
#include "radioactive.h"

using std::cout;
using std::cerr;
using std::endl;
namespace pc = physical_constants;

//-----------------------------------------------------------------
// set opacity and emissivities
//-----------------------------------------------------------------
void transport::set_opacity(double dt)
{

  double tend,tstr;
  double get_system_time(void);

  // tmp vector to hold emissivity
  vector<OpacityType> emis(nu_grid_.size());
  vector<OpacityType> scat(nu_grid_.size());
  vector<OpacityType> res(nu_grid_.size());
  vector<OpacityType> fluor(nu_grid_.size());
  vector<OpacityType> bb_emis(nu_grid_.size());
  emis.assign(emis.size(),0.0);

  // always do LTE on first step

  if (first_step_) {
    for (auto i_gas_state = gas_state_vec_.begin(); i_gas_state != gas_state_vec_.end(); i_gas_state++) {
      if (first_step_) i_gas_state->use_nlte_ = 0;
    }
  }

  // zero out opacities, etc...
  for (int i=0;i<grid->n_zones;i++)
  {
    compton_opac[i]  = 0;
    photoion_opac[i] = 0;
    rosseland_mean_opacity_[i] = 0;
    planck_mean_opacity_[i]    = 0;
    emissivity_[i].wipe();
    for (int j=0;j<nu_grid_.size();j++)
    {
      abs_opacity_[i][j] = 0;
      if (!omit_scattering_) scat_opacity_[i][j] = 0;
      if (use_line_res_scatter_) {
        line_res_opacity_[i][j] = 0;
        line_fluor_opacity_[i][j] = 0;
      }
    }
  }


  if (verbose)
    if (solve_Tgas_with_updated_opacities_ && first_step_ == 0)
      printf("# Solving coupled equations for gas state and temperature\n");


  // loop over my zones to calculate
  // loop to parallelize with OpenMP
  tstr = get_system_time();
  int solve_root_errors = 0;
  int solve_iter_errors = 0;
  int solve_Tbracket_errors = 0;
  int solve_Titer_errors = 0;

#pragma omp parallel firstprivate(emis, scat, res, fluor, bb_emis) shared(cerr,solve_root_errors,solve_iter_errors,solve_Tbracket_errors,solve_Titer_errors) default(none)
  {
#ifdef _OPENMP
    int my_threadID = omp_get_thread_num();
#else
    int my_threadID = 0;
#endif
    // Private variables for each thread
    GasState* gas_state_ptr = &(gas_state_vec_[my_threadID]);
    vector<double> X_now(grid->n_elems);
    radioactive radio_obj;
    radioactive* radio = &radio_obj;
    int solve_error = 0;


#pragma omp for reduction(+:solve_Tbracket_errors,solve_Titer_errors)
    for (int i=my_zone_start_;i<my_zone_stop_;i++) {
      // pointer to current zone for easy access
      zone* z = &(grid->z[i]);

      //------------------------------------------------------
      // calculate optical photon opacities
      //------------------------------------------------------

      // set up the state of the gas in this zone
      gas_state_ptr->dens_ = z->rho;
      gas_state_ptr->temp_ = z->T_gas;
      gas_state_ptr->time_ = t_now_;
      if (gas_state_ptr->temp_ < temp_min_value_) gas_state_ptr->temp_ = temp_min_value_;
      if (gas_state_ptr->temp_ > temp_max_value_) gas_state_ptr->temp_ = temp_max_value_;

      // Compute local |dv/dr| for Sobolev opacity (non-monotonic velocity support).
      // The standard Sobolev approximation uses 1/t_exp as the velocity gradient,
      // which is wrong for non-monotonic velocity (e.g., shock interaction).
      // Here we compute the actual local dv/dr from the grid velocity profile.
      // For 1D spherical grid, dv/dr = (v[i] - v[i-1]) / (r[i] - r[i-1]).
      // At velocity extrema (dv/dr ~ 0), fall back to 1/t_exp to avoid
      // infinite Sobolev depth.
      {
        double dvds_sob = 0;
        double v_this = grid->z[i].v[0];
        // innermost zone: velocity is extrapolated flat inside rmin
        // (v_inner_ = z[0].v[0]), so dv/dr = 0 -> fallback to 1/t_exp
        double v_prev = (i > 0) ? grid->z[i-1].v[0] : v_this;
        double dr;
        grid->get_zone_size(i, &dr);
        if (dr > 0) {
          dvds_sob = fabs(v_this - v_prev) / dr;
        }
        // If dv/dr is too small (near velocity extremum), use 1/t_exp
        double dvds_homologous = 1.0 / t_now_;
        if (dvds_sob < 0.1 * dvds_homologous) {
          dvds_sob = dvds_homologous;
        }
        gas_state_ptr->sobolev_dvds_ = dvds_sob;
      }

      // radioactive decay the composition
      for (size_t j=0;j<X_now.size();j++) X_now[j] = z->X_gas[j];
      if (!omit_composition_decay_) {
        radio->decay_composition(grid->elems_Z,grid->elems_A,X_now,t_now_);
      }

      gas_state_ptr->set_mass_fractions(X_now);

      gas_state_ptr->bulk_grey_opacity_ = z->bulk_grey_opacity;
      gas_state_ptr->total_grey_opacity_ = z->total_grey_opacity;

      // NLTE density threshold: in low-density zones the rate matrix becomes
      // singular (collisional rates ~ne -> 0, radiative rates ~J_nu -> 0,
      // all multiplied by lev_lte_ floored at 1e-30). Fall back to LTE for
      // zones below nlte_dens_min_ to avoid GSL LU decomposition crash.
      int nlte_was_on = gas_state_ptr->use_nlte_;
      if (nlte_was_on && z->rho < nlte_dens_min_)
        gas_state_ptr->use_nlte_ = 0;

      // Lazy allocation of the LTE solve cache (one entry per zone)
      if (lte_cache_valid_.size() != (size_t)grid->n_zones)
      {
        lte_state_cache_.resize(grid->n_zones);
        lte_cache_T_.assign(grid->n_zones, 0.0);
        lte_cache_rho_.assign(grid->n_zones, 0.0);
        lte_cache_X_.assign(grid->n_zones, vector<double>());
        lte_cache_valid_.assign(grid->n_zones, 0);
      }

      if (first_step_)
      {
        zone* z = &(grid->z[i]);
        gas_state_ptr->dens_ = z->rho;
        gas_state_ptr->temp_ = z->T_gas;

        if (gas_state_ptr->total_grey_opacity_ == 0)
        {
          // always do LTE on first step, without updating temperature
          solve_error = gas_state_ptr->solve_state();
        }
      }

      else
      {
        if (solve_Tgas_with_updated_opacities_)
        {
          // gas state solution (LTE or NLTE) solution as well as radiative
          // equilibrium temperature solve will happen here
          solve_error = solve_state_and_temperature(gas_state_ptr, i);
        }

        else
        {
          if (gas_state_ptr->total_grey_opacity_ == 0)
          {
            // LTE solves are deterministic in (rho, T, X) and ignore
            // J_nu, so a cached state with identical inputs is exact.
            bool cache_hit = false;
            if (!gas_state_ptr->use_nlte_ && lte_cache_valid_[i]
                && z->rho == lte_cache_rho_[i]
                && gas_state_ptr->temp_ == lte_cache_T_[i]
                && X_now == lte_cache_X_[i])
            {
              *gas_state_ptr = lte_state_cache_[i];
              cache_hit = true;
            }
            if (!cache_hit)
            {
              solve_error = gas_state_ptr->solve_state(J_nu_[i]);
              if (!gas_state_ptr->use_nlte_)
              {
                lte_state_cache_[i] = *gas_state_ptr;
                lte_cache_T_[i]   = gas_state_ptr->temp_;
                lte_cache_rho_[i] = z->rho;
                lte_cache_X_[i]   = X_now;
                lte_cache_valid_[i] = 1;
              }
            }
          }
        }
      }

      if (solve_error == 1) solve_root_errors += 1;
      if (solve_error == 2) solve_iter_errors += 1;
      if (solve_error == 3) solve_Tbracket_errors += 1;
      if (solve_error == 4) solve_Titer_errors += 1;

      //gas_state_ptr->print();
      #pragma omp critical
      if(write_levels) gas_state_ptr->write_levels(i);

      grid->z[i].n_elec = gas_state_ptr->n_elec_;

      // restore NLTE flag if it was temporarily disabled for low-density zone
      gas_state_ptr->use_nlte_ = nlte_was_on;

      // calculate the opacities/emissivities
      gas_state_ptr->computeOpacity(abs_opacity_[i],scat,emis,res,fluor,bb_emis);

      double max_extinction = maximum_opacity_* z->rho;

      // save and normalize emissivity cdf
      grid->z[i].L_thermal = 0;
      if (nu_grid_.size() == 1)
      {
        double bb_int = pc::sb*pow(grid->z[i].T_gas,4)/pc::pi;
        grid->z[i].L_thermal += 4*pc::pi*abs_opacity_[i][0]*bb_int;
        emissivity_[i].set_value(0,1);
        if (!omit_scattering_) scat_opacity_[i][0] = scat[0];
      }
      else for (int j=0;j<nu_grid_.size();j++)
      {
        double ednu = emis[j]*nu_grid_.delta(j);
        emissivity_[i].set_value(j,ednu);
        grid->z[i].L_thermal += 4*pc::pi * ednu;
        if (!omit_scattering_) scat_opacity_[i][j] = scat[j];
        if (use_line_res_scatter_) {
          line_res_opacity_[i][j] = res[j];
          line_fluor_opacity_[i][j] = fluor[j];
          line_emis_cdf_[i].set_value(j, bb_emis[j]*nu_grid_.delta(j));
        }

        // check for maximum opacity
        if (!omit_scattering_)
        {
          if (scat_opacity_[i][j] > max_extinction)
            scat_opacity_[i][j]   = max_extinction;
        }
        if (abs_opacity_[i][j]  > max_extinction)
          abs_opacity_[i][j]    = max_extinction;
        if (use_line_res_scatter_)
        {
          if (line_res_opacity_[i][j] > max_extinction)
            line_res_opacity_[i][j] = max_extinction;
          if (line_fluor_opacity_[i][j] > max_extinction)
            line_fluor_opacity_[i][j] = max_extinction;
        }
      }
      emissivity_[i].normalize();
      if (use_line_res_scatter_) line_emis_cdf_[i].normalize();

      // flattened importance CDF: sample_j = (1-f)*emis_j + f/nbins,
      // weight_j = emis_j/sample_j. Emitted energy is unchanged in the
      // mean; rare line bins just get sampled far more often.
      if (emit_flatten_ > 0)
      {
        int nb = nu_grid_.size();
        double f = emit_flatten_;
        for (int j=0;j<nb;j++)
        {
          double e = emissivity_[i].get_value(j);
          double s = (1.0 - f)*e + f/nb;
          emit_cdf_[i].set_value(j,s);
        }
        emit_cdf_[i].normalize();
        for (int j=0;j<nb;j++)
        {
          double s = emit_cdf_[i].get_value(j);
          emit_wt_[i][j] = (s > 0) ? emissivity_[i].get_value(j)/s : 0;
        }
      }

      // Optional per-zone emissivity/opacity dump for diagnostics.
      // SEDONA_DUMP_EMIS_ZONES="i1 i2 ..." writes one text file per
      // zone per opacity call: nu, emissivity density, abs, scat.
      if (getenv("SEDONA_DUMP_EMIS_ZONES"))
      {
        static std::map<int,int> emis_call;
        std::istringstream iss(getenv("SEDONA_DUMP_EMIS_ZONES"));
        int zid;
        while (iss >> zid)
        {
          if (zid != i) continue;
          int call = ++emis_call[i];
          char fname[128];
          sprintf(fname,"emis_dump_z%04d_c%02d.txt",i,call);
          FILE *ef = fopen(fname,"w");
          fprintf(ef,"# zone %d call %d  T=%.4g rho=%.3g  L_th=%.3e\n",
              i,call,grid->z[i].T_gas,z->rho,grid->z[i].L_thermal);
          for (int j=0;j<nu_grid_.size();j++)
            fprintf(ef,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                nu_grid_.center(j), emis[j],
                abs_opacity_[i][j], scat_opacity_[i][j],
                use_line_res_scatter_ ? line_res_opacity_[i][j] : 0,
                use_line_res_scatter_ ? line_fluor_opacity_[i][j] : 0);
          fclose(ef);
        }
      }

      // calculate mean opacities. The resonance channel is real
      // extinction for diffusion, so it enters the means together
      // with the thermalizing part of the line opacity.
      vector<OpacityType> abs_for_means;
      if (use_line_res_scatter_)
      {
        abs_for_means = abs_opacity_[i];
        for (int j=0;j<nu_grid_.size();j++)
          abs_for_means[j] += line_res_opacity_[i][j] + line_fluor_opacity_[i][j];
      }
      const vector<OpacityType>& abs_mean =
        use_line_res_scatter_ ? abs_for_means : abs_opacity_[i];
      planck_mean_opacity_[i] =
        gas_state_ptr->get_planck_mean(abs_mean,scat_opacity_[i]);
      rosseland_mean_opacity_[i] =
        gas_state_ptr->get_rosseland_mean(abs_mean,scat_opacity_[i]);
      if (ddmc_n_groups_ > 1)
        gas_state_ptr->get_group_rosseland_means(abs_mean,
            scat_opacity_[i], nu_to_group_, rosseland_group_opacity_[i]);

      //------------------------------------------------------
      // gamma-ray opacity (compton + photo-electric)
      //------------------------------------------------------
      compton_opac[i]  = 0;
      photoion_opac[i] = 0;
      for (int k=0;k<grid->n_elems;k++)
      {
        double dens  = z->X_gas[k]*z->rho;
        double ndens = dens/(pc::m_p*grid->elems_A[k]);
        // compton scattering opacity
        compton_opac[i] += ndens*pc::thomson_cs*grid->elems_Z[k];
        // photoelectric opacity
        double photo = pow(pc::alpha_fs,4.0)*4.0*sqrt(2.0);
        photo *= pow(1.0*grid->elems_Z[k],5.0);
        photo *= pow(pc::m_e_MeV,3.5);
        photoion_opac[i] += ndens*2.0*pc::thomson_cs*photo;
      }
    }

    // output any solve error
    #pragma omp single
    if (verbose) 
    {
      if (solve_root_errors != 0) 
        std::cerr << "# Warning: root not bracketed in n_e solve in " << solve_root_errors << " zones" << std::endl;
      if (solve_iter_errors != 0) 
        std::cerr << "# Warning: max iterations hit in n_e solve in " << solve_iter_errors << " zones" << std::endl;
      if (solve_Tbracket_errors != 0)
        std::cerr << "# Warning: root not bracketed in T solve in " << solve_Tbracket_errors << " zones (kept previous T)" << std::endl;
      if (solve_Titer_errors != 0)
        std::cerr << "# Warning: max iterations hit in T solve in " << solve_Titer_errors << " zones (kept previous T)" << std::endl;
    }
  }
  // end OpenMP parallel region
 
  if (solve_Tgas_with_updated_opacities_ && first_step_ == 0) {
    reduce_Tgas(); }


  reduce_n_elec();

  tend = get_system_time();
  if (verbose) cout << "# Calculated opacities   (" << (tend-tstr) << " secs) \n";


  //------------------------------------------------------------
  // Calcuate implicit MC parameter eps_imc
  //------------------------------------------------------------
  for (int i=0;i<grid->n_zones;i++)
  {
    if (radiative_eq)
    {
      grid->z[i].eps_imc = 1.;
    }
    else
    {
      // Not distinguishing between lab frame density and comoving frame density
      double fleck_beta  = 4.0*pc::a*pow(grid->z[i].T_gas,4)/(grid->z[i].e_gas*grid->z[i].rho);
      // here planck mean opac has units cm^-1 .
      // When grey opacity is used, planck_mean_opacity should just be the correct grey opacity
      double tfac = pc::c*planck_mean_opacity_[i]*dt;
      double f_imc = fleck_alpha_*fleck_beta*tfac;

      // make sure to avoid divide by zero if fleck alpha is zero and not computing e_gas through hydro
      if (fleck_alpha_ == 0)
        f_imc = 0.;

      grid->z[i].eps_imc = 1.0/(1.0 + f_imc);
    }
  }

  // turn nlte back on after first step, if wanted
  if (first_step_) {
    for (auto i_gas_state = gas_state_vec_.begin(); i_gas_state != gas_state_vec_.end(); i_gas_state++) {
      if (first_step_) i_gas_state->use_nlte_ = use_nlte_;
    }
  }
}


//-----------------------------------------------------------------
// get comoving opacity at the frequency
// returns the frequency index of the photon in the
// comoving frame
//-----------------------------------------------------------------
int transport::get_opacity(particle &p, double dshift, double &opac, double &eps)
{
  double res_frac, fluor_frac;
  return get_opacity(p,dshift,opac,eps,res_frac,fluor_frac);
}

int transport::get_opacity(particle &p, double dshift, double &opac, double &eps,
                           double &res_frac, double &fluor_frac)
{
  assert(p.ind >= 0);
  res_frac = 0;
  fluor_frac = 0;

  // comoving frame frequency
  double nu = p.nu*dshift;
  int i_nu = 0;

  // get opacity if it is an optical photon.
  if (p.type == photon)
  {
    // interpolate opacity at the local comving frame frequency
    i_nu = nu_grid_.locate_within_bounds(nu);
    double a_opac = nu_grid_.value_at(nu,abs_opacity_[p.ind],i_nu);
    double s_opac = 0;
    if (!omit_scattering_) s_opac = nu_grid_.value_at(nu,scat_opacity_[p.ind],i_nu);
    double l_opac = 0, f_opac = 0;
    if (use_line_res_scatter_)
    {
      l_opac = nu_grid_.value_at(nu,line_res_opacity_[p.ind],i_nu);
      f_opac = nu_grid_.value_at(nu,line_fluor_opacity_[p.ind],i_nu);
    }
    opac = a_opac + s_opac + l_opac + f_opac;
    if (opac == 0) eps = 0;
    else
    {
      eps        = a_opac/opac;
      res_frac   = l_opac/opac;
      fluor_frac = f_opac/opac;
    }
  }

  // get opacity if it is a gamma-ray
  if (p.type == gammaray)
  {
    double c_opac = compton_opac[p.ind]*klein_nishina(p.nu);
    double p_opac = photoion_opac[p.ind]*pow(p.nu,-3.5);
    opac = c_opac + p_opac;
    eps  = p_opac/(c_opac + p_opac);
  }

  return i_nu;
}


//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double transport::klein_nishina(double x)
{
  // divide by m_e c^2 = 0.511 MeV
  x = x/pc::m_e_MeV;
  double logfac = log(1 + 2*x);
  double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
  double term2 = 1.0/2.0/x*logfac;
  double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
  double KN    = .75*(term1 + term2 + term3);
  return KN;
}


//-----------------------------------------------------------------
// calculate planck function in frequency units
//-----------------------------------------------------------------
double transport::blackbody_nu(double T, double nu)
{
  double zeta = pc::h*nu/pc::k/T;
  return 2.0*nu*nu*nu*pc::h/pc::c/pc::c/(exp(zeta)-1);
}
