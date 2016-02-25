

//  Created by Heather Ratcliffe on 11/02/2016.


#ifndef _non_thermal_h
#define _non_thermal_h


class non_thermal{

private:
  bool configure_from_file(std::string file_prefix);

public:

  non_thermal(std::string file_prefix);
  ~non_thermal();
  calc_type ref_dens;
  calc_type ref_B;
  calc_type fraction;
  calc_type v_par;
  calc_type v_perp;


};





#endif