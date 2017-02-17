

//  Created by Heather Ratcliffe on 11/02/2016.


#ifndef _non_thermal_h
#define _non_thermal_h
#include<map>
#include "support.h"
/** \brief Nonthermal electron description
*
*Small class to hold a non-thermal electron distribution we can operate on and with. Distribution is defined by reading a <filepath>nonthermal.conf file. This either specifies a functional form and corresponding constants as references into the deck.status file, or a lookup file containing a data_array. The interface remains the same, accessing either f_p(p_par, p_perp) or df/dp.
*
*An example conf file is \verbatim
ncomps = 1
nonrely = 1
hot:
function = max
dens=dens_rath
vpar =vtherm_parh
vperp=vtherm_perph
end:
\endverbatim specifying a single Bimaxwellian with density dens_rath from deck.status, etc. Note the "end:" line. Any number of components with any form can be supplied and the resulting f_p is their sum. The nonrely flag requests to use the non-relativistic calculations. An example of a lookup based nonthermal.conf is \verbatim
ncomps = 1
hot:
function =lookup
lookup = my_data.dat
end: 
\endverbatim where the file my_data contains a data_array. For example, the output of compress_distributions utility can be used, or the IDL routines in refit_distribs.
*
*To add new functional forms, create a function such as bimax, below, and create the binding in non_thermal::configure_from_file under "Binding function free parameters to create f_p"
*/

class non_thermal{

private:
  std::vector<std::function<calc_type(calc_type p_par, calc_type p_perp)> > f_p_private;/**< Bound function(s) for f*/
  std::map<std::string, calc_type> parameters;/**<Store for deck.status file parameters*/
  calc_type total_dens;/**<Total summed density of all components */
  calc_type dp;/**<delta p for calculating df/dp*/
  bool norely;/**<Flag for non-relativistic calculation (changes normalisation of Maxwellians)*/

/********Basic setup and allocation functions ****/
  bool configure_from_file(std::string file_prefix);

public:

  my_type * lookup_data;/** Data pointer for use with a lookup type function backend. Note type matched to MY EPOCH data*/
  size_t dims[2];

/********Basic setup and allocation functions ****/
  non_thermal(std::string file_prefix);
  ~non_thermal();/**<Clean up. Calls clean_lookup() which can be used to do anything needed to cleanup after a lookup function*/

/********Primary interface functions ****/
  calc_type f_p(calc_type p_par, calc_type p_perp);
  calc_type d_f_p(calc_type p_par, calc_type p_perp, bool parallel);

/********Helper functions ****/
  void set_dp(calc_type dp){this->dp = dp;}/**<Set the dp used to get numerical derivative*/
  calc_type get_total_dens(){return this->total_dens;}/**<Return the total density of components*/
  bool get_norely(){return this->norely;}/**<Return flag for non-relativistic calculation*/

/********Special testing functions ****/
#ifdef RUN_TESTS_AND_EXIT
  calc_type get_v_par(){if(parameters.count("vtherm_par")>0)return parameters["vtherm_par"]; else if(parameters.count("vtherm_parh")>0)return parameters["vtherm_parh"]; else return 0.0;}
  calc_type get_v_perp(){if(parameters.count("vtherm_perp")>0)return parameters["vtherm_perp"]; else if(parameters.count("vtherm_perph")>0)return parameters["vtherm_perph"]; else return 0.0;}
  calc_type get_t_par(){/**WORKS FOR SINGLE COMP MAX ONLY*/ if(parameters.count("temp_par")>0)return parameters["temp_par"];else
    return 0.0;}
  calc_type get_t_perp(){/**WORKS FOR SINGLE COMP MAX ONLY*/ if(parameters.count("temp_perp")>0)return parameters["temp_perp"];else
  return 0.0;}
#endif

};

/********Functions providing component specifications ****/
calc_type bikappa(calc_type p, calc_type p2, calc_type kappa, calc_type v_k, calc_type v_k2, calc_type A);

calc_type bimax(calc_type p, calc_type p2, calc_type p_th, calc_type p_th2, calc_type A);

calc_type lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min);
calc_type seperable_lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min);

/* *******Helpers for component functions ****/

std::function<calc_type(calc_type p_par, calc_type p_perp)> configure_lookup(std::string file_prefix,std::string file, non_thermal * my_nonth);
void clean_lookup(non_thermal & my_nonth);

#endif