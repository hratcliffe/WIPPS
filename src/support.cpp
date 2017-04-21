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
#include <regex>

#include "support.h"
#include "data_array.h"
#include "reader.h"

mpi_info_struc mpi_info = mpi_info_null;/**< Global mpi_info structure*/

deck_constants my_const;/**< Global deck_constants*/

//----------- SPECIFIC HELPER FUNCTIONS -------------

/********MPI and code helpers ****/
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
  
  
  MPI_Type_free(&deckConstType);

}

/********IO helpers ****/
void get_deck_constants(std::string file_prefix){
/** \brief Setup run specific constants
*
*Reads deck.status and parses values for user defined constants etc. It will rely on using the specific deck, because it has to look for names. Any changes to deck may need updating here. Tag names are set as const strings in support.h. IMPORTANT: If we find additional density tags we fold those into om_pe. To prevent this, either remove from deck.status, or prefix their printed names with something so they do not match the strings in support.h
@param file_prefix File prefix prepended to "deck.status"
*/


//open file deck.status from whatever path
//Find the line " Constant block values after"
//Read into vector up to "Deck state" or eof
//parse out the values we want given by name in support.h

  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string line;
  std::vector<std::string> lines;
  bool found=false;
  
  my_const.omega_ce = 0.0;
  my_const.omega_pe = 0.0;
  my_const.dens_factor = 0.0;
  my_const.omega_ci = 0.0;
  my_const.v_t = 0.0;
  my_const.ppc = 0;
  
  if(!infile.is_open()){
    my_error_print("No deck.status file found, aborting", mpi_info.rank);
    return;
  }
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS.size()).c_str(), CONSTANTS.c_str()))){
      found = true;
      break;
    }
  }
  if(!found){
    my_error_print("No constants dump found in deck.status, aborting", mpi_info.rank);
    return;
  }

  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS_END.size()).c_str(), CONSTANTS_END.c_str()))) break;
    lines.push_back(line);
  }
  
  std::string name, val;
  bool parse_err;
  float val_f;
  
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
    my_const.dens_factor = 1.0+ tmp_rat+tmp_rath;
    
  }
  my_const.omega_ci = my_const.omega_ce * me/mp;
  //assumes same charge magnitude, I.e. H plasma

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
    if(strcmp(argv[i], "-f")==0 && i < argc-1){
      values.file_prefix = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-start")==0 && i < argc-1){
      if(atoi(argv[i+1]) >= 0){
        values.time[0] = atoi(argv[i+1]);
      }else{
        my_error_print("Start cannot be negative!!!!!!");
      }
      i++;
    }
    else if(strcmp(argv[i], "-end")==0 && i < argc-1){
      if(atoi(argv[i+1]) >= 0){
        values.time[1] = atoi(argv[i+1]);
      }else{
        my_error_print("End cannot be negative!!!!!");
      }
      i++;
    }
    else if(strcmp(argv[i], "-rows")==0 && i < argc-1){
      if(atoi(argv[i+1]) >= 0){
        values.time[2] = atoi(argv[i+1]);
      }else{
        my_error_print("Rows cannot be negative!!!!!");
      }
      i++;
    }
    else if(strcmp(argv[i], "-block")==0 && i < argc-1){
      values.block = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-n")==0 && i < argc-1){
      if(atoi(argv[i+1]) >= 0){
        values.n_space = atoi(argv[i+1]);
      }else{
        my_error_print("N cannot be negative!!!!!");
      }
      i++;
    }
    else if(strcmp(argv[i], "-space")==0 && i < argc-2){
      if(atoi(argv[i+1]) >= 0 && atoi(argv[i+2]) >= 0){
        values.space[0] = atoi(argv[i+1]);
        values.space[1] = atoi(argv[i+2]);
      }else{
        my_error_print("Space cannot be negative!!!!!");
      }
      i+=2;
    }
    else if(strcmp(argv[i], "-d")==0 && i < argc-2){
      if(atoi(argv[i+1]) >= 0 && atoi(argv[i+2]) >= 0){
        values.d[0] = atoi(argv[i+1]);
        values.d[1] = atoi(argv[i+2]);
      }else{
        my_error_print("D cannot be negative!!!!!");
      }
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
    else if(!((strlen(argv[i]) > 0) && argv[i][0] == HANDLED_ARG[0])){
      my_print("UNKNOWN OPTION "+mk_str(argv[i]), 0);
    }
  }

  if(values.space[0] == -1 && values.space[1] == -1 && values.n_space == -1){
    values.n_space = mpi_info.n_procs;
    values.per_proc = 1;
    //By default do one block per processor
  }
  else if(values.n_space == -1){
    //If we've asked for just a single block with given limits, do that. We will NOT account for asking for one block on more than 1 processor
    values.n_space = 1;
    values.per_proc = 1;
  }
  else{
    values.per_proc = values.n_space / mpi_info.n_procs;
    //integer division!
    values.n_space = values.per_proc * mpi_info.n_procs;
    //Exactly divide
  }

  if(values.time[1] < values.time[0]) values.time[1] = values.time[0] + 1;
  //If unspecified use 1 row per file
  if(values.d[0] >MAX_SIZE){
    values.d[0] = MAX_SIZE;
    my_error_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  if(values.d[1] >MAX_SIZE){
    values.d[1] = MAX_SIZE;
    my_error_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  //Protect from invalid user input
  
  return values;
}

void process_command_line_help_arg(int argc, char *argv[], char help_id){
/** \brief Handle help request at command line
*
* If -h is supplied anywhere in arg list, print the appropriate help file and then exit
*/

  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-h")==0){
      print_help(help_id);
      exit(0);
    }
  }
}

void print_help(char code){
/** \brief Print command line help
*
*Prints contents of halp_file from rank zero and calls safe exit.
@param code Input single character utility name code to get specific help. If non-empty file opened is halp_file with '_'+code inserted before extension
*/
  std::ifstream halp;
  
  std::string file = halp_file;
  if(code !=0 && code !=' ') file = append_into_string(file, "_"+std::string(1, code));
  halp.open(file);
  if(!halp){
    //File not found!
    my_print("Help file "+file+" not found", 0);
    return;
  }
  if(mpi_info.rank == 0 && halp){
    std::cout<<halp.rdbuf();
    std::cout<<'\n';
  }
}

void log_code_constants(std::string file_prefix){
/** \brief Log internal constants
*
*Records ID codes etc as name value pairs
@param file_prefix File path to write to
\todo Updates?
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
    file<<"DEFAULT_SPECTRUM_ANG_STDDEV "<<DEFAULT_SPECTRUM_ANG_STDDEV<<'\n';
    file<<"SPECTRUM_THRESHOLD "<<SPECTRUM_THRESHOLD<<'\n';

  }
  file.close();
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

/********Data helpers ****/
void divide_domain(std::vector<size_t> dims, size_t space[2], int per_proc, int block_num){
/** \brief Divide dims evenly between procs
*
*Uses the number of space blocks from args (if specified) and the domain size from dims to ensure perfect subdivision and set current proc's bounds. We can ignore incoming space vals as they should be -1
@param dims Vector of sizes to divide
@param[out] space 2-element array giving the local space section for this processor
@param per_proc Number of blocks per processor See setup_args per_proc field or see process_command_line() for example of calculation
@param block_num Current block index on this processor
*/

    size_t end, block_start, block_end, block_len;
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

my_type get_ref_Bx(std::string file_prefix, size_t space_in[2], size_t time_0, bool is_acc){
/** Read reference B_x from the specfied file prefix as given, dump number time_0
@param file_prefix File path
@param space_in Limits on x-dimension to slice out
@param time_0 The dump time to read
@param is_acc Whether these files use the accumulation extension to EPOCH
@return Average bx over specified space range at given time
*/
  data_array bx = get_Bx(file_prefix, space_in, time_0, is_acc);
  if(bx.is_good()) return avval(bx);
  else return 0.0;
}

data_array get_Bx(std::string file_prefix, size_t space_in[2], size_t time_0, bool is_acc){
/** Read reference B_x from file at path file_prefix, dump number time_0. If space_in is not [-1, -1], only the slice it dictates is read
@param file_prefix File path
@param space_in Limits on x-dimension to slice out
@param time_0 The dump time to read
@param is_acc Whether these files use the accumulation extension to EPOCH
@return data_array containing bx data
*/
  my_print("Getting ref B");
  char block_id[ID_SIZE];
  if(!is_acc) strcpy(block_id, "bx");
  else strcpy(block_id, "abx");
  
  reader bx_reader = reader(file_prefix, block_id);
  //We use this to get the local average B field
  size_t bx_times[3] = {time_0, time_0+1, 1};
  //use specified file and read one row
  bx_reader.update_ref_filenum(time_0);

  size_t n_dims;
  std::vector<size_t> dims;
  int err = bx_reader.read_dims(n_dims, dims);
  if(err) return data_array();

  size_t space_dim = space_in[1]-space_in[0];
  data_array bx = data_array(space_dim, 1);
  if(!bx.is_good()) return data_array();
  
  if(n_dims == 1){
    err = bx_reader.read_data(bx, bx_times, space_in);
  }else if(n_dims == 2){
    err = bx_reader.read_data(bx, bx_times, space_in, 1);
  }else{
    my_error_print("3-D space not added...", mpi_info.rank);
  }
  if(err == 1 || err == 2) return bx;
  //2 is a non-fatal read error
  else return data_array();
}

bool flatten_fortran_slice(my_type * src_ptr, my_type* dest_ptr, size_t n_dims_in, size_t * dims_in, size_t flatten_on_dim, size_t flat_start, size_t flat_stop){
/** \brief Flatten a Fortran-style array on the specified dimension
*
* The result is a Fortran-style array of rank n_dims_in - 1, containing the total along each value of the flattening dim. dest_ptr is assumed to point to an allocated block sufficient to hold the result. NB this produces a total not an average. NB flat_start and flat_stop must be valid indices into the dim to be flattened. They CANNOT be checked. If supplied only the slice they delimit is totalled
@param src_ptr Pointer to start of source data
@param dest_ptr Pointer to start of destination memory block
@param n_dims_in Rank (number of dimensions) of input "array"
@param dims_in Dimensions of input
@param flatten_on_dim Dimension to flatten on
@param flat_start Start of slice to flatten (default 0)
@param flat_stop End of slice to flatten (default MAX_SIZE_T)
@return 0 for success, 1 if parameters invalid (usually a range error)
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
  
  size_t els_to_copy = 1, n_segments = 1, sz = 1;

  for(size_t i=0; i<flatten_on_dim; ++i) els_to_copy*= dims_in[i];
  for(size_t i=flatten_on_dim+1; i< n_dims_in; ++i) n_segments *= dims_in[i];

  size_t chunk_sz = els_to_copy;
  for(size_t i=0; i< n_segments; ++i) std::copy(src_ptr + i*chunk_sz*dims_in[flatten_on_dim] +els_to_copy*flat_start, src_ptr + i*chunk_sz*dims_in[flatten_on_dim]+ els_to_copy+ els_to_copy*flat_start, dest_ptr + i*chunk_sz*sz);
  
  //Now we should have the 0th row in place
  //Add each successive row onto it
  
  for(size_t j = flat_start+1; j<flat_stop; j++){
    
    for(size_t i=0; i< n_segments; ++i) std::transform(src_ptr + i*chunk_sz*dims_in[flatten_on_dim] + els_to_copy*j, src_ptr + i*chunk_sz*dims_in[flatten_on_dim]+ els_to_copy + els_to_copy*j, dest_ptr + i*chunk_sz*sz,dest_ptr + i*chunk_sz*sz, std::plus<my_type>());
  }

  return 0;
  
}

int where(const my_type * ax_ptr, int len, my_type target){
/** \brief Find where ax_ptr exceeds target
*
*Checks bounds and calls locally scoped recursive whereb to do the search. Special case if target equals bottom end
@param ax_ptr Pointer to start of axis to search
@param len Length of axis we're searching
@param target Target value to find
@return Index of target in axis
*/
  int whereb(const my_type * ax_ptr, int len, my_type target, int &cut,int sign); //Recursive function to do the finding

  int sign = 1.0;
  int cut = 0;
  if(ax_ptr[0] == target) return 0;
  if(ax_ptr[0] > target || ax_ptr[len-1] < target) return -1;
  return whereb(ax_ptr, len, target, cut, sign);

}

int whereb(const my_type * ax_ptr, int len, my_type target,int &cut, int sign){
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

//----------- HELPER TYPE FUNCTION DECLARATIONS -------------

/********IO helpers ****/
std::string read_wipps_version_string(std::string filename){
/** \brief Read code version from file
*
*Reads the version string of code used to write the file given by (full path) filename
@param filename File to read
@return String containing code version
*/

  std::fstream infile;
  infile.open(filename, std::ios::in|std::ios::binary);
  if(!infile.good()){
    my_error_print("File open or access error");
    return "";
  }
  char tmp_vers[GIT_VERSION_SIZE];
  //Skip past the 3 numbers before string
  infile.seekg(2*sizeof(int)+sizeof(my_type), std::ios::beg);
  infile.read((char*) &tmp_vers, sizeof(char)*GIT_VERSION_SIZE);
  infile.close();
  return std::string(tmp_vers);
}

bool check_wipps_version(std::string filename){
/** \brief Check wipps versioning
*
*Reads version from specified file and compares, printing error and returning result.
@param filename File to read
@return 1 if match, 0 else
*/

  std::string wipps_version = read_wipps_version_string(filename);
  if(!compare_as_version_string(wipps_version)){
    my_error_print("Warning, a different code version was used to write this file. Data may be incompatible");
    return 0;
  }
  return 1;
}

void my_print(std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
*MPI aware screen output. Prints from one or all processors to cout
@param text Text to print
@param rank Rank of this processor
@param rank_to_write Which rank should do the printing, default 0. Set to -1 to print from all
@param noreturn Set to not output a line break after text
*/
  if(rank == rank_to_write || rank_to_write == -1){
  
    std::cout<< text;
    if(!noreturn)std::cout<<std::endl;
  }

}
void my_print(std::fstream * handle, std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
*MPI aware file output. Prints from one or all processors to given filestream
@param handle Filestream to write to. If this is nullptr, cout is used
@param text Text to print
@param rank Rank of this processor
@param rank_to_write Which rank should do the printing, default 0. Set to -1 to print from all
@param noreturn Set to not output a line break after text
*/
  if((rank == rank_to_write || rank_to_write == -1) && handle!=nullptr){
    *handle<<text;
    if(!noreturn)*handle<<std::endl;
  }else if(rank == rank_to_write || rank_to_write == -1){
    std::cout<<text;
    if(!noreturn)std::cout<<std::endl;

  }

}

void my_error_print(std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
*MPI aware screen error output. Prints from one or all processors to cerr
@param text Text to print
@param rank Rank of this processor
@param rank_to_write Which rank should do the printing, default 0. Set to -1 to print from all
@param noreturn Set to not output a line break after text
*/
  if(rank == rank_to_write || rank_to_write == -1){
  
    std::cerr<< text;
    if(!noreturn)std::cerr<<std::endl;
  }

}
void my_error_print(std::fstream * handle, std::string text, int rank, int rank_to_write, bool noreturn){
/** \brief Write output
*
*MPI aware filestream error output. Prints from one or all processors to given filestream. If this is nullptr, cerr is used
@param handle Filestream to write to. If this is nullptr, cout is used
@param text Text to print
@param rank Rank of this processor
@param rank_to_write Which rank should do the printing, default 0. Set to -1 to print from all
@param noreturn Set to not output a line break after text
*/
  if((rank == rank_to_write || rank_to_write == -1) && handle!=nullptr){
    *handle<<text;
    if(!noreturn)*handle<<std::endl;
  }else if(rank == rank_to_write || rank_to_write == -1){
    std::cerr<<text;
    if(!noreturn)std::cerr<<std::endl;

  }

}

/********String handling helpers ****/
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
std::string mk_str(long double i, bool noexp){return mk_str((double) i, noexp);}
std::string mk_str(float i, bool noexp){return mk_str((double) i, noexp);}
std::string mk_str(char * str){
  return std::string(str);
}
void trim_string(std::string &str, char ch){
/** \brief Trim ch from ends of string
*
*Trims leading and trailing occurences of character from string
@param str The string (will be changed)
@param ch The character to trim
*/
  std::string tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(str.find_first_not_of(ch), str.size());
  str=tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(0, str.find_last_not_of(ch)+1);
  str=tmp;

}

std::string replace_char(std::string str_in, char ch, char repl){
/** \brief Replace character in string
*
*Replace one character with another in a string
@param str_in Input string
@param ch Character to replace
@param repl Character to replace with
@return The amended string
*/

  std::string str = str_in;
  size_t pos =str.find_first_of(ch);
  while(pos != std::string::npos){
    str[pos] = repl;
    pos =str.find_first_of(ch);
  }
  return str;
}

std::string append_into_string(const std::string &in, const std::string &infix){
/** \brief Insert infix in string
*
*Inserts the infix string into in BEFORE the last file extension. If no '.' is found in string, append to end. First char being . is not an extension.
@param in String to insert into. Usually ends in .[extension]
@param infix String to insert
@return The resulting modified string
*/
  size_t start = 0;
  if(in.size() > 0) start = in.substr(1, in.size()).find_last_of('.') +1;

  std::string in_copy = in;
  if(start !=std::string::npos){
    in_copy.insert(start, infix);
    return in_copy;
  }else{
    return in+infix;
  }
}

bool parse_name_val(std::string in, std::string &name, std::string &val){
/** \brief Parse x=y strings
*
* Basic line parser. Takes a string and if it contains an '=' splits into the left and right segments, stripping leading and trailing spaces. Returns 0 if success, 1 if no equals sign. Standard comment character is # as first non-whitespace
@param in The input string to parse
@param[out] name The name part
@param[out] val The value part
@return 0 if line parsed, 1 if it can't be
*/

  if(in.find_first_not_of(" \t\n") == std::string::npos) return 1;
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

bool compare_as_version_string(std::string str, std::string vers_str, bool minor){
/** \brief Check str against version code
*
*Checks string against git version code VERSION. Version strings are vx.y if present at all. By default check just the x, if minor is true, check y also. If either string doesn't match expected format we try comparing as just strings. This maintains behaviour from before I used version tags where just commit-id was checked
@param str String to check
@param vers_str String to check against, defaults to VERSION
@param minor Flag to check minor version number too
@return True if equal, false else
*/

  std::string major_v="0", minor_v="0", major_in="0", minor_in="0";
  std::regex re("v([0-9]+).([0-9]+)");
  std::smatch matches;
  bool bad_version=false, bad_input=false;
  if(std::regex_search(vers_str, matches, re) && matches.size() > 2){
    major_v = matches[1];
    minor_v = matches[2];
  }else{
    bad_version = true;
  }
  if(std::regex_search(str, matches, re) && matches.size() > 2){
    major_in = matches[1];
    minor_in = matches[2];
  }else{
    bad_input = true;
  }
  if(bad_version || bad_input){
    //Try comparing as just strings.
    return str == vers_str;
  }
  if(!minor){
    return major_v == major_in;
  }else{
    return (major_v == major_in) && (minor_v == minor_in);
  }
}

/********Maths helpers ****/
template<typename T> T integrator(T * start, int len, T * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory.
@param start Pointer to start of data
@param len Length of data
@param increment Pointer to start of increment axis (e.g. x[1:end]-x[0:end-1])
@return The integrated value
*/

  T value = 0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }

 return value;

}

template float integrator<float>(float *, int, float *);
template double integrator<double>(double *, int, double *);
//We need both float and double versions

calc_type square_integrator(calc_type * start, int len, calc_type * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory.
@param start Pointer to start of data
@param len Length of data
@param increment Pointer to start of increment axis (e.g. x[1:end]-x[0:end-1])
@return The square-integrated value
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
@param start Start of data interval to smooth
@param len Length of interval to smooth
@param width Boxcar smoothing width
@param periodic Flag to set for wrap-around smoothing
*/

  if(width > len) my_print("Smoothing width exceeds array size", mpi_info.rank);
  calc_type result = 0.0;
  int edge = width/2;
  int off = width - edge*2;
  T * hold;
  hold = (T*) malloc(len*sizeof(T));
  for(int i=0; i<width; ++i) result += start[i];
  //Sum first width elements
  for(int i=edge+off; i<len-edge; ++i){
    result -= start[i-edge-off];
    result += start[i+edge];
    //remove behind and add in front. Faster for width>2
    hold[i] = result / (calc_type) width;
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
  //Copy back into place
  std::copy(hold, hold+len, start);
  free(hold);

}
template void inplace_boxcar_smooth(float *, int, int, bool);
template void inplace_boxcar_smooth(double *, int, int, bool);

std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn){
/** \brief Finds roots of cubic x^3 + an x^2 + bn x + cn = 0
*
* Uses Num. Rec. equations, which are optimised for precision. Note that if x >>1 precision errors may result. Returns real solutions only
@param an Coefficient of x^2
@param bn Coefficient of x
@param cn Constant part
@return Vector of real roots
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
  return ret_vec;

}

//[2] below tells compiler nothing but is helpful to remember
template<typename T> T interpolate_linear(T axis[2], T vals[2], T target){
/** Interpolate vals on axis to target value. 
*
*Axis and vals should contain 2 values boxing the target. We use linear interpolation to obtain the axis value corresponding to target 
@param axis Axis values for interpolation 
@param vals Values at axis values 
@param target Target value to interpolate to 
@return The interpolated value
*/

  T ret = (std::abs(target - axis[1]) * vals[0] + std::abs(target - axis[0]) * vals[1])/(std::abs(axis[1] - axis[0]));
  return ret;
}

template float interpolate_linear(float*, float*, float);
template double interpolate_linear(double*, double*, double);
//Again we need both float and double versions

//[2] below tells compiler nothing but is helpful to remember
template<typename T> T interpolate_nearest(T axis[2], T vals[2], T target){
/** Interpolate vals on axis to target value. 
*
*Axis and vals should contain 2 values boxing the target. We select the nearest axis value as corresponding to target 
@param axis Axis values for interpolation 
@param vals Values at axis values 
@param target Target value to interpolate to 
@return The interpolated value
*/

  T ret = 0.0;

  //select closer value
  if(std::abs(target - axis[0]) <= std::abs(target - axis[1])) ret = vals[0];
  else ret = vals[1];
 
  return ret;
}
template float interpolate_nearest(float*, float*, float);
template double interpolate_nearest(double*, double*, double);
//Again we need both float and double versions
