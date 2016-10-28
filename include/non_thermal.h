

//  Created by Heather Ratcliffe on 11/02/2016.


#ifndef _non_thermal_h
#define _non_thermal_h
#include<map>
#include "support.h"
/** \brief Nonthermal electron description
*
*Very small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status assuming particular input deck format. We use a std::function which allows us to bind any number of constants to values and then invoke the function using single interface. That is f_p_par(p) can map to exp(-p^2/const) or equally (1+1/const*p^kappa)^1/kappa. Alternately a lookup table for numerically defined function can be used behind-the-scenes*/

class non_thermal{

private:
  bool configure_from_file(std::string file_prefix);
  std::vector<std::function<calc_type(calc_type p_par, calc_type p_perp)> > f_p_private;/**< Bound function for f*/
  std::map<std::string, calc_type> parameters;
  std::map<std::string, std::string> extra_parameters;
  calc_type norm;
  size_t ncomps;
  calc_type total_dens;
  calc_type dp;
  calc_type ref_dens;/**<Reference density (background)*/
  calc_type ref_B;/**<Reference B field*/
  calc_type fraction;/**<Non-thermal fraction*/
  calc_type v_par;/**<Parallel velocity*/
  calc_type v_perp;/**Perpendicular velocity \todo Clean up redundant params*/
  bool norely;
public:

  non_thermal(std::string file_prefix);
  ~non_thermal(){if(lookup_data) free(lookup_data);};
  void write(std::ofstream &outfile);
  void dump(std::fstream &outfile);
  calc_type f_p(calc_type p_par, calc_type p_perp);
  calc_type d_f_p(calc_type p_par, calc_type p_perp, bool parallel);
  size_t get_n_pars(){return parameters.size();}
  void set_dp(calc_type dp){this->dp = dp;}/**<Set the dp used to get numerical derivative*/
  my_type * lookup_data;/** Data pointer for use with a lookup type function backend. Note type matched to MY EPOCH data*/
//#ifdef RUN_TESTS_AND_EXIT
  calc_type get_v_par(){return v_par;}
  calc_type get_v_perp(){return v_perp;}
//#endif
  calc_type get_t_par(){/**WORKS FOR SINGLE COMP MAX ONLY*/ if(parameters.count("temp_par")>0)return parameters["temp_par"];else
    return v_par*v_par*me/kb;}
  calc_type get_t_perp(){/**WORKS FOR SINGLE COMP MAX ONLY*/ if(parameters.count("temp_perp")>0)return parameters["temp_perp"];else
  return v_perp*v_perp*me/kb;}
  calc_type get_total_dens(){return this->total_dens;}
  bool get_norely(){return this->norely;}
};

calc_type bikappa(calc_type p, calc_type p2, calc_type kappa, calc_type v_k, calc_type v_k2, calc_type A);
calc_type max(calc_type p, calc_type v_th, calc_type A);
calc_type double_max(calc_type p, calc_type p_th, calc_type p_th2, calc_type A, calc_type A2);
calc_type bimax(calc_type p, calc_type p2, calc_type p_th, calc_type p_th2, calc_type A);

calc_type lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min);
calc_type seperable_lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min);

std::function<calc_type(calc_type p_par, calc_type p_perp)> configure_lookup(std::string file_prefix,std::string file, non_thermal * my_nonth);



#endif