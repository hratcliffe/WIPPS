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


class spectrum;
class diffusion_coeff;

/** \brief Plasma parameters and dispersion
*
*This class will take care of solving the plasma dispersion roots, selecting a wave mode if necessary etc etc. The base calculations are modified from file mufunctions3.f90 author  Clare E. J. Watt                                   date 18/05/10. From there: Note that mu is calculated using the Appleton-Hartree relation, and the choice of sign is obtained from Albert [2005].
* @author Heather Ratcliffe  @date 07/10/2015
*/
class plasma{

private:
  static const int ncomps=4; 

  calc_type B0;/**<Local (block averged) B0*/
  calc_type om_ce_local;/**<Local actual cyclotron frequency*/
  calc_type om_ce_ref;/**< REFERENCE cyclotron freq*/

  //Other Plasma params, such as species massm charge
  calc_type pmass[ncomps];
  calc_type pcharge[ncomps];
  calc_type pdens[ncomps];

  bool configure_from_file(std::string file_prefix);
  bool is_setup;
public:

  plasma(){is_setup = false;}
  plasma(std::string file_prefix, my_type Bx_local=-1);
  ~plasma();

  bool is_good(){return is_setup;}/**<Whether everything is setup*/
  void get_density(){;}/**< Density is assumed constant*/
  void get_B0(){;}/**<B0 can vary in space*/

  void set_B0(my_type B0);
  
  mu get_root(calc_type th, calc_type w, calc_type psi, bool Righthand=true);
  mu_dmudom get_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n=0, bool Righthand=true)const;
  mu_dmudom get_high_dens_phi_mu_om(calc_type w, calc_type psi, calc_type alpha, int n, calc_type omega_n, bool Righthand=true)const;

  std::vector<calc_type> get_resonant_omega(calc_type x, calc_type v_par, calc_type n)const;

  calc_type get_omega_ref(std::string code)const;
  calc_type get_dispersion(my_type k, int wave_type, bool reverse=0, bool deriv=0, my_type theta=0.0)const;

};

#endif
