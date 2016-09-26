//
//  main_support.cpp
//  
//
//  Created by Heather Ratcliffe on 25/02/2016.
//
//

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <map>
#include <cstdlib>

#include "support.h"
#include "main_support.h"
#include "data_array.h"
#include "reader.h"

mpi_info_struc mpi_info = mpi_info_null;/**< MPI data Not const as defined from output of MPI_Init. Use "extern const mpi_info_struc mpi_info;
" to access elsewhere. This may be terrible, but this doesn't warrant a class, really, come on.*/

extern deck_constants my_const;

/** \defgroup main Main Helper Functions
*@{ */

int local_MPI_setup(int argc, char *argv[]){
/** \brief Do the MPI init
*
* Calls MPI init and sets up communicator. Stores number of processors and each rank into mpi_info.
*/
  int ierr, rank, n_procs;

  ierr = MPI_Init(&argc, &argv);
  //Note any other command line arg processing should account for this...

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  mpi_info.rank = rank;
  mpi_info.n_procs = n_procs;
  
  //if(rank ==0) std::cout<< &mpi_info<<std::endl;

  
  return ierr;
}

void safe_exit(){
/** \brief Exit program
*
* Does minimal cleanup and exits
*/

  if(mpi_info.n_procs >0) MPI_Finalize();

  exit(0);
}

void share_consts(){
/** \brief MPI Share deck constants
*
*Share the deck constants read on root to all procs
*/

  const int count = 5;
  int lens[count] ={1, 1, 1, 1, 1};
  MPI_Aint disps[count];
  MPI_Datatype types[count] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
  MPI_Datatype deckConstType;
  
  disps[0] = offsetof(deck_constants, v_t);
  disps[1] = offsetof(deck_constants, omega_pe);
  disps[2] = offsetof(deck_constants, omega_ce);
  disps[3] = offsetof(deck_constants, omega_ci);
  disps[4] = offsetof(deck_constants, ppc);
  
  MPI_Type_create_struct(count, lens, disps, types, &deckConstType);
  MPI_Type_commit(&deckConstType);

  MPI_Bcast(&my_const, 1, deckConstType, 0, MPI_COMM_WORLD);

}

std::vector<std::string> process_filelist(int argc, char *argv[]){
/** \brief Extracts files from list
*
*  Returns vector of filename strings. Assumed to be in form [stuff]_space0_space1.dat and will be ordered on space0
*/

  std::vector<std::string> names;
  for(int i=0; i< argc; i++){
    if(((strcmp(argv[i], "-Finput")==0)||(strcmp(argv[i], "-Sinput")==0)) && i < argc-1){
      while(i<argc-1 && argv[i+1][0]!= '-'){
        //Checks if next argument is a new flag
        names.push_back(argv[i+1]);
        i++;
      }
    }
  }

  //Now we have vector of names. Extract space0 from them and plop into a map
  std::map<int, std::string> name_map;
  int num;
  size_t posa, posb;
  for(size_t i = 0; i< names.size(); i++){
    posb = names[i].find_last_of('_');
    posa = (names[i].substr(0, posb)).find_last_of('_');
    num = atoi((names[i].substr(posa+1, posb-posa-1)).c_str());
    name_map[num] = names[i];
  }
  names.clear();
  //Now we plop the map back into the vector
  for(auto it=name_map.begin(); it!=name_map.end(); it++){
    names.push_back(it->second);
  }
  //This ends up with sorting without having to write custom comparator. Replace when have time
  
  return names;
}

setup_args process_command_line(int argc, char *argv[]){
/** \brief Set basic parameters
*
*Sets defaults or those given via command line (see help.txt). Forces an constant integer number of space blocks on each core.
*/

  setup_args values;
  values.n_space = -1;
  values.space[0] = -1;
  values.space[1] = -1;

  values.time[0] = 0;
  values.time[1] = 1;
  values.time[2] = 0;
  values.use_row_time = false;
  values.file_prefix = "./files/";
  values.block = "ex";
  values.d[0] = 10;
  values.d[1] = 10;
  values.is_list = false;
  values.is_spect = false;

  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-h")==0){
      print_help();
      exit(0);
    }
    else if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-start")==0 && i < argc-1){
      values.time[0] = atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-end")==0 && i < argc-1){
      values.time[1] = atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-rows")==0 && i < argc-1){
      values.time[2] = atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-block")==0 && i < argc-1){
      values.block = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-n")==0 && i < argc-1){
      values.n_space= atoi(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-space")==0 && i < argc-2){
      values.space[0] = atoi(argv[i+1]);
      values.space[1] = atoi(argv[i+2]);
      i+=2;
    }
    else if(strcmp(argv[i], "-d")==0 && i < argc-2){
      values.d[0] = atoi(argv[i+1]);
      values.d[1] = atoi(argv[i+2]);
      i+=2;
    }
    else if(((strcmp(argv[i], "-Finput")==0)||(strcmp(argv[i], "-Sinput")==0)) && i < argc-1){
      values.is_list=true;
      //Now hunt for next arg..., we assume no '-' starting filenames
      int tmp = i;
      while(i<argc-1 && argv[i+1][0]!= '-') i++;
      if(tmp -i >= 1 ) i--;
      //Go back one so that loop advance leaves us in correct place, but not if we didn't skip on at all or we'd infinite loop
    }
    else if((strcmp(argv[i], "-Sinput")==0) && i < argc-1){
      values.is_spect = true;
      //Now hunt for next arg..., we assume no '-' starting filenames
      int tmp = i;
      while(i<argc-1 && argv[i+1][0]!= '-') i++;
      if(tmp -i >= 1 ) i--;
    }
    else if(!(strlen(argv[i]) > 0) && argv[i][0] ==HANDLED_ARG[0]) std::cout<<"UNKNOWN OPTION " <<argv[i]<<'\n';

  }

  if(values.space[0] == -1 && values.space[1] == -1 && values.n_space == -1){
    values.n_space = mpi_info.n_procs;
    values.per_proc = 1;
  }
  //By default do one block per processor
  else{
    values.per_proc = values.n_space / mpi_info.n_procs;
    //integer division!
    values.n_space = values.per_proc * mpi_info.n_procs;
  //Exactly divide
  }

  if(values.time[0]< 0 ) values.time[0] = 0;
  if(values.time[1]< 0 ) values.time[1] = 0;
  if(values.time[1]< 0 ) values.time[1] = 0;
  if(values.time[1] < values.time[0]) values.time[1] = values.time[0] + 1;
  //If unspecified use 1 row per file
  if(values.d[0] < 0) values.d[0] = 0;
  if(values.d[1] < 0) values.d[1] = 0;
  if(values.d[0] >MAX_SIZE){
    values.d[0] = MAX_SIZE;
    my_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  if(values.d[1] >MAX_SIZE){
    values.d[1] = MAX_SIZE;
    my_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  //Protect from invalid user input
  
  return values;
}

void print_help(char code){
/** \brief Print command line help
*
*Prints contents of halp_file from rank zero and calls safe exit. Input single character utility name code to get specific help (assumed to be in halp_file_[c]
*/
  std::ifstream halp;
  
  std::string file = halp_file;
  if(code !=0) file = append_into_string(file, "_"+std::string(1, code));
  halp.open(file);
  if(!halp){
    //File not found!
    std::cout<<"Help file "<<file<<" not found\n";
    return;
  }
  if(mpi_info.rank == 0 && halp){
    std::cout<<halp.rdbuf();
    std::cout<<'\n';
  }
}

void divide_domain(std::vector<size_t> dims, int space[2], int per_proc, int block_num){
/** \brief Divide dims evenly between procs
*
*Uses the number of space blocks from args (if specified) and the domain size from dims to ensure perfect subdivision and set current proc's bounds. We can ignore incoming space vals as they should be -1
*/

    int end, block_start, block_end, block_len;
    float per_proc_size;
    end = dims[0];
    per_proc_size = std::ceil( (float) end / (float) mpi_info.n_procs);
    //Force overlap not missing
    block_start = mpi_info.rank * (int) std::floor((float) end/(float) (mpi_info.n_procs));
    block_end = block_start + (int) per_proc_size;
    block_len = block_end - block_start;
    space[0] = block_start + (block_num) * block_len / per_proc;
    space[1] = space[0] + block_len / per_proc;

}

int where(my_type * ax_ptr, int len, my_type target){
/** \brief Find where ax_ptr exceeds target
*
*Checks bounds and calls locally scoped recursive whereb to do the search. Special case if target equals bottom end
*/
  int whereb(my_type * ax_ptr, int len, my_type target, int &cut,int sign); //Recursive function to do the finding

  int sign = 1.0;
  int cut = 0;
  if(ax_ptr[0] == target) return 0;
  if(ax_ptr[0] > target || ax_ptr[len-1] < target) return -1;
  return whereb(ax_ptr, len, target, cut, sign);

}

int whereb(my_type * ax_ptr, int len, my_type target,int &cut, int sign){
/**\brief Recursive binary bisection find. 
*
*First index where ax_ptr exceeds target. If sign ==-1 it should do < but as yet untested....
*/
  
  if(len==1){
     cut++;
     return cut;
  }
  else if(sign*ax_ptr[len/2] == sign*target){
    //Exact equality. Save recursing deeper
    cut+= len/2;
    return cut;
  }
  else if(sign*ax_ptr[len/2] > sign*target){
    whereb(ax_ptr, len/2, target, cut, sign);

    return cut;
  }
  else if(sign*ax_ptr[len/2] < sign*target){
    cut+= len/2;
    whereb(ax_ptr+len/2, len-len/2, target, cut, sign);
    return cut;
  }

  return -1;
  //Shouldn't ever reach this case, but.
}

void get_deck_constants(std::string file_prefix){
/** \brief Setup run specific constants
*
*Reads deck.status and parses values for user defined constants etc. It will rely on using the specific deck, because it has to look for names. Any changes to deck may need updating here. Tag names are set as const strings in support.h. IMPORTANT: If we find additional density tags we fold those into om_pe. To prevent this, either remove from deck.status, or prefix their printed names with something so they do not match the strings in support.h
*/


//open file deck.status from whatever path
//Find the line " Constant block values after"
//Read into vector up to "Deck state" or eof
//parse out the values we want given by name in support.h

  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string header_row, line;
  infile>> header_row;
  std::vector<std::string> lines;
  bool found=false;
  if(!infile.is_open()){
    my_print("No deck.status file found, aborting", mpi_info.rank);
    safe_exit();
  
  }
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS.size()).c_str(), CONSTANTS.c_str()))){
      found = true;
      break;
    }
  }
  if(!found){
    my_print("No constants dump found in deck.status, aborting", mpi_info.rank);
    safe_exit();
  }

  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS_END.size()).c_str(), CONSTANTS_END.c_str()))) break;
    lines.push_back(line);
  }
  
  std::string name, val;
  bool parse_err;
  float val_f;
  
  my_const.omega_ce = 0;
  my_const.omega_pe = 0;
  float tmp_rat=0.0, tmp_rath=0.0;
  
  for(size_t i=0; i< lines.size(); i++){

    parse_err = parse_name_val(lines[i], name, val);
    if(parse_err) continue;
    val_f = atof(val.c_str());
  
    if(name == OMEGA_CE) my_const.omega_ce = val_f;
    else if(name == OMEGA_PE) my_const.omega_pe = val_f;
    else if(name == DENS_RAT) tmp_rat = val_f;
    else if(name == DENS_RATH) tmp_rath = val_f;
  }

  if(tmp_rat != 0.0 || tmp_rath !=0.0){
    my_print("Modifying density to " + mk_str(1.0+tmp_rat+tmp_rath, true));
    my_const.omega_pe *= std::sqrt(1.0+tmp_rat+tmp_rath);
    
  }
  my_const.omega_ci = my_const.omega_ce * me/mp;
  //assumes same charge magnitude, I.e. H plasma

}

/** @} */

void log_code_constants(std::string file_prefix){
/** \brief Log internal constants
*
*Records ID codes etc as name value pairs \todo Updates?
*/
  std::ofstream file;
  std::string filename = file_prefix +"constants.dump";
  file.open(filename.c_str());
  if(file.is_open()){
  
    file<<"WAVE_WHISTLER "<<WAVE_WHISTLER<<'\n';
    file<<"WAVE_PLASMA "<<WAVE_PLASMA<<'\n';
    file<<"WAVE_O "<<WAVE_O<<'\n';
    file<<"FUNCTION_NULL "<<FUNCTION_NULL<<'\n';
    file<<"FUNCTION_DELTA "<<FUNCTION_DELTA<<'\n';
    file<<"FUNCTION_GAUSS "<<FUNCTION_GAUSS<<'\n';
    file<<"FUNCTION_ISO "<<FUNCTION_ISO<<'\n';
    file<<"V_MIN "<<V_MIN<<'\n';
    file<<"V_MAX "<<V_MAX<<'\n';
    file<<"ANG_MIN "<<ANG_MIN<<'\n';
    file<<"ANG_MAX "<<ANG_MAX<<'\n';
    file<<"SPECTRUM_ANG_STDDEV "<<SPECTRUM_ANG_STDDEV<<'\n';
    file<<"SPECTRUM_THRESHOLD "<<SPECTRUM_THRESHOLD<<'\n';

  }
  file.close();
}

bool parse_name_val(std::string in, std::string &name, std::string &val){
/** \brief Parse x=y strings
*
* Basic line parser. Takes a string and if it contains an '=' splits into the left and right segments, stripping leading and trailing spaces. Returns 0 if success, 1 if no equals sign. Standard comment character is # as first non-whitespace
*/
  if(in[in.find_first_not_of(" \t\n")] == '#') return 1;
  size_t pos = in.find("=");
  if(pos == std::string::npos){
    //is not a x = y line
    name = "";
    val="";
    return 1;

  }else{
    name = in.substr(0, pos);
    trim_string(name);
    val = in.substr(pos+1, in.size());
    trim_string(val);
    return 0;
  }
}

void trim_string(std::string &str, char ch){
  std::string tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(str.find_first_not_of(ch), str.size());
  str=tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(0, str.find_last_not_of(ch)+1);
  str=tmp;

}

std::string replace_char(std::string str_in, char ch, char repl){
  std::string str = str_in;
  size_t pos =str.find_first_of(ch);
  while(pos != std::string::npos){
    str[pos] = repl;
    pos =str.find_first_of(ch);
  }
  return str;
}
void my_print(std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also. Set rank_to_write to -1 to dump from all. Default value is 0
*/
  if(rank == rank_to_write || rank_to_write == -1){
  
    std::cout<< text;
    if(!noreturn)std::cout<<std::endl;
  }

}
void my_print(std::fstream * handle, std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also. Set rank_to_write to -1 to dump from all. Default value is 0
*/
  if((rank == rank_to_write || rank_to_write == -1) && handle!=nullptr){
    *handle<<text;
    if(!noreturn)*handle<<std::endl;
  }else if(rank == rank_to_write || rank_to_write == -1){
    std::cout<<text;
    if(!noreturn)std::cout<<std::endl;

  }

}

std::string append_into_string(const std::string &in, const std::string &infix){
/** \brief Insert infix in string
*
*Inserts the infix string into in BEFORE the last file extension. If no '.' is found in string, append to end. First char being . is not an extension.
*/
  size_t start = in.substr(1, in.size()).find_last_of('.') +1;
  std::string in_copy = in;
  if(start !=std::string::npos){
    in_copy.insert(start, infix);
    return in_copy;
  }else{
    return in+infix;
  }
}

std::string mk_str(int i){

  char buffer[25];
  std::sprintf(buffer, "%i", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(size_t i){

  char buffer[25];
  std::sprintf(buffer, "%lu", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(double i, bool noexp){

  char buffer[25];
  if(noexp) std::snprintf(buffer, 25, "%f", i);
  else std::snprintf(buffer, 25, "%e", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(bool b){

  if(b) return "1";
  else return "0";

}

std::string mk_str(long double i, bool noexp){return mk_str((double) i, noexp);};
std::string mk_str(float i, bool noexp){return mk_str((double) i, noexp);};

template<typename T> T integrator(T * start, int len, T * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  T value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template float integrator<float>(float *, int, float *);
template double integrator<double>(double *, int, double *);
//We need both float and double versions

calc_type square_integrator(calc_type * start, int len, calc_type * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  calc_type value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i]*start[i] + start[i+1]*start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template<typename T> void inplace_boxcar_smooth(T * start, int len, int width, bool periodic){
/** \brief Boxcar smoothing of specified width
*
*Smooths the array given by start and len using specified width. If periodic is set the ends wrap around. Otherwise they one-side
*/

  if(width > len) my_print("Really? How am I meant to smooth that?", mpi_info.rank);
  calc_type result = 0.0;
  int edge = width/2;
  for(int i=0; i<width; ++i) result += start[i];
  for(int i=edge+1; i<len-edge; ++i){
    result -= start[i-edge-1];
    result += start[i+edge];
    //remove behind and add in front. Faster for width>2
    start[i] = result / (calc_type) width;
  }

  //Handle ends
  if(periodic){
    result = 0.0;
    for(int i=0; i<width-1; ++i) result += start[i];
    //first width-1
    int wrap = len -1;
    //and one wrapped back
    result +=start[wrap];
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - i - 1];
      start[edge-i] = result / (calc_type) width;
    }
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - edge - i - 1];
      start[len-i-1] = result / (calc_type) width;
    }


  }else{
    //we just truncate at the bottom
    
  
  }


}
template void inplace_boxcar_smooth(calc_type *, int, int, bool);
template void inplace_boxcar_smooth(my_type *, int, int, bool);

std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn){
/** \brief Finds roots of cubic x^3 + an x^2 + bn x + cn = 0
*
* Uses Num. Rec. equations, which are optimised for precision. Note that if x >>1 precision errors may result. Returns real solutions only
*/

  calc_type Q, R, bigA, bigB, Q3, R2, bigTheta;
  std::vector<calc_type> ret_vec;

  Q = (std::pow(an, 2) - 3.0 * bn)/9.0;
  R = (2.0* std::pow(an, 3) - 9.0 * an *bn + 27.0*cn)/54.0;
  
  R2 = std::pow(R, 2);
  Q3 = std::pow(Q, 3);
  
  if( R2 < Q3){
    
    bigTheta = std::acos(R/sqrt(Q3));
    calc_type minus2sqrtQ = -2.0*std::sqrt(Q);
    
    ret_vec.push_back(minus2sqrtQ*std::cos(bigTheta/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta + 2.0*pi)/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta - 2.0*pi)/3.0) - an/3.0);

  }else{
    calc_type ret_root;
    bigA = - boost::math::sign(R)*std::pow((std::abs(R) + std::sqrt(R2 - Q3)), 1.0/3.0 );

    (bigA != 0.0) ? (bigB = Q / bigA) : (bigB = 0.0);
    ret_root = (bigA + bigB) - an/3.0;

    ret_vec.push_back(ret_root);
  }

/** Used to test when writing
  calc_type tmp;
  for(int i=0; i<ret_vec.size(); ++i){
    
    tmp = std::pow(ret_vec[i], 3) + an*std::pow(ret_vec[i], 2) + bn*ret_vec[i] + cn;
    std::cout<<"solution gives "<<tmp<<std::endl;
  
  }
*/
  return ret_vec;

}

template<typename T> T interpolate(T* axis, T* vals, T target, int pts){
/** Interpolate vals on axis to target value
*
*For pts=1 uses closest value, pts=2 uses 2 pt linear, \todo add more pts options
*/

  T ret = 0.0;
  if(pts ==1){
    //select closer value
    if(std::abs(target - axis[0]) <= std::abs(target - axis[1])) ret = vals[0];
    else ret = vals[1];
  
  }else if(pts ==2){
    ret = (std::abs(target - axis[1]) * vals[0] + std::abs(target - axis[0]) * vals[1])/(std::abs(axis[1] - axis[0]));
  
  }

  return ret;
}
template float interpolate(float*, float*, float, int);
template double interpolate(double*, double*, double, int);

my_type get_ref_Bx(std::string file_prefix, int space_in[2], int time_0, bool is_acc){
/** Read reference B_x from the specfied file prefix as given, dump number time_0*/
  my_print("Getting ref B");
  char block_id[ID_SIZE];
  if(!is_acc) strcpy(block_id, "bx");
  else strcpy(block_id, "abx");
  
  reader bx_reader = reader(file_prefix, block_id);
  //We use this to get the local average B field
  int bx_times[3] = {time_0, time_0+1, 1};
  //use specified file and read one row
  bx_reader.set_ref_filenum(time_0);

  size_t n_dims;
  std::vector<size_t> dims;
  int err = bx_reader.read_dims(n_dims, dims);
  if(err) return 0.0;

  size_t space_dim = space_in[1]-space_in[0];
  data_array bx = data_array(space_dim, 1);
  if(!bx.is_good()) return 0.0;
  
  if(n_dims == 1){
    err = bx_reader.read_data(bx, bx_times, space_in);
  }else if(n_dims == 2){
    err = bx_reader.read_data(bx, bx_times, space_in, 1);
  }else{
    my_print("3-D space not added...", mpi_info.rank);
  }
  
  if(err == 0 || err ==2 ) return bx.avval();
  //2 is a non-fatal read error
  else return 0.0;
}

bool flatten_fortran_slice(my_type * src_ptr, my_type* dest_ptr, size_t n_dims_in, size_t * dims_in, size_t flatten_on_dim, size_t flat_start, size_t flat_stop){
/** \brief Flatten a Fortran-style array on the specified dimension
*
* The result is a Fortran-style array of rank n_dims_in - 1, containing the total along each value of the flattening dim. dest_ptr is assumed to point to an allocated block sufficient to hold the result. NB this produces a total not an average
*/
/*
* A 2-d array 5x3 is
|ooooo||ooooo||ooooo|
<-row->

A 3-d 5x3x2 is
[|ooooo||ooooo||ooooo|][|ooooo||ooooo||ooooo|]
<--------'slice'------>

We'll borrow the resizer code to get the new array by resizing the required dim to size 1 in a copy, and then we'll add the rest on
Think of the array as being 3-D. The dim we;re flattening is dim-1. All the less-significant dims are smushed together into dim-0, and all the higher into dim-2
*/

  if(flatten_on_dim > n_dims_in) return 1;
  if(flat_start >= dims_in[flatten_on_dim]) return 1;
  if(flat_stop >= dims_in[flatten_on_dim]) flat_stop = dims_in[flatten_on_dim];
  
  //my_print("Flattening array", mpi_info.rank);

  size_t part_sz = 1;
  size_t els_to_copy = 1, n_segments = 1, sz = 1;

  for(size_t i=0; i<flatten_on_dim; ++i) els_to_copy*= dims_in[i];
  for(size_t i=flatten_on_dim+1; i< n_dims_in; ++i) n_segments *= dims_in[i];

  part_sz = els_to_copy*n_segments;

  size_t chunk_sz = els_to_copy;
  for(size_t i=0; i< n_segments; ++i) std::copy(src_ptr + i*chunk_sz*dims_in[flatten_on_dim] +els_to_copy*flat_start, src_ptr + i*chunk_sz*dims_in[flatten_on_dim]+ els_to_copy+ els_to_copy*flat_start, dest_ptr + i*chunk_sz*sz);
  
  //Now we should have the 0th row in place
  //Add each successive row onto it
  
  for(size_t j = flat_start+1; j<flat_stop; j++){
    
    for(size_t i=0; i< n_segments; ++i) std::transform(src_ptr + i*chunk_sz*dims_in[flatten_on_dim] + els_to_copy*j, src_ptr + i*chunk_sz*dims_in[flatten_on_dim]+ els_to_copy + els_to_copy*j, dest_ptr + i*chunk_sz*sz,dest_ptr + i*chunk_sz*sz, std::plus<my_type>());
  }

  return 0;
  
}