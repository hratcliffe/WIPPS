//
//  plasma.h
//  
//
//  Created by Heather Ratcliffe on 07/10/2015.
//
//

#ifndef _plasma_h
#define _plasma_h

#include<vector>
#include "support.h"
#include "resonance_poly.h"


class spectrum;
class diffusion_coeff;

typedef enum p_state {p_none, p_good, p_default, p_overflow, p_underflow} plasma_state;/**< State of plasma and result of attempted define from file*/

/** \brief Plasma parameters and dispersion
*
*Plasma objects contain specifications for a plasma, including density, B field, and species composition. In general we configure them from a file, plasma.conf, and no other constructor is provided. After construction, the plasma will be valid, but if given file is not found or reading fails, default values will be used. Two cyclotron frequencies are available, local and reference. For varying B fields, a reference B field should be given and the local om_ce will match this. The reference value always matches that in the deck_constants struct.
*We offer functions to solve plasma dispersion "exactly" using the various get_root, get_phi* etc functions, or get_dispersion which uses various analytic approximations, usually high density ones. For the high-density approximations, the FIRST species is assumed to be the electrons. The former functions are modified from file mufunctions3.f90 author  Clare E. J. Watt (18/05/10). From there: Note that mu is calculated using the Appleton-Hartree relation, and the choice of sign is obtained from Albert\cite Albert2005.
*An additional function to simultaneously solve the Doppler type resonance condition, and the approximate dispersion relation are provided.
*IMPORTANT: Plasma setup relies on my_consts being defined!!
* @author Heather Ratcliffe  @date 07/10/2015 \ingroup cls
* \caveat Note that no spatial variations in density or ambient B field are included in the dispersion calculations here

*/

class plasma{

private:
  static const int ncomps=4;/**<Number of component species*/

  calc_type B0;/**<Local (block averged) B0*/
  calc_type om_ce_local;/**<Local actual cyclotron frequency*/
  calc_type om_ce_ref;/**< REFERENCE cyclotron freq*/

  calc_type pmass[ncomps];/**<Per species mass*/
  calc_type pcharge[ncomps];/**<Per species charge*/
  calc_type pdens[ncomps];/**<Per species density*/
  calc_type pvth[ncomps];/**<Per species thermal velocity*/

  //resonance_poly * full_poly;
  NR_poly * full_poly;

  plasma_state is_setup;/**<Check for validity*/
/********Basic setup and allocation functions ****/
  plasma_state configure_from_file(std::string file_prefix);

/********Dispersion solver core ****/  
  mu_dmudom get_phi_mu_from_stix(calc_type w, calc_type psi,calc_type alpha, int n, calc_type gamma_particle, calc_type R, calc_type L, calc_type P, bool skip_phi, bool Righthand)const;

public:

/********Basic setup and allocation functions ****/
  explicit plasma(){is_setup = p_none;}/**<Default constructor, create useless plasma object*/
  explicit plasma(std::string file_prefix, my_type Bx_local=-1);

/********Get/set functions ****/
  bool is_good(){return is_setup != p_none;}/**<Whether everything is setup @return Boolean true is good, false else*/
  calc_type get_omega_ref(std::string code)const;
  calc_type get_B0(){return B0;}/**<Return B0. This can vary in space @return Value of B0 for this plasma*/
  void set_B0(my_type B0);

/********Dispersion solvers ****/
  /** \brief Solve plasma dispersion only
  *
  * Solve dispersion, omitting extended phi calcs
  @param w Wave frequency
  @param psi Wave normal angle 
  @return mu_dmudom struct containing mu and derivs */
  mu_dmudom get_mu(calc_type w, calc_type psi) const{return get_phi_mu_om(w, psi, 0.0, 0, 1.0, true);}
  /** \brief Solve plasma dispersion only
  *
  * Solve dispersion, omitting extended phi calcs, using high_dens approximation
  @param w Wave frequency
  @param psi Wave normal angle 
  @return mu_dmudom struct containing mu and derivs */
  mu_dmudom get_high_dens_mu(calc_type w, calc_type psi) const{return get_high_dens_phi_mu_om(w, psi, 0.0, 0, 1.0, true);}
  mu_dmudom get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, bool skip_phi = false, bool Righthand=true)const;
  mu_dmudom get_high_dens_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type gamma_particle, bool skip_phi=false, bool Righthand=true)const;

  std::vector<calc_type> get_resonant_omega(calc_type theta, calc_type v_par, calc_type gamma_particle, int n)const;
  std::vector<calc_type> get_resonant_omega_full(calc_type theta, calc_type v_par, calc_type gamma_particle, int n)const;
  calc_type get_dispersion(my_type k, int wave_type, bool reverse=0, bool deriv=0, my_type theta=0.0)const;
// calc_type gamma_particle,
};

#endif
