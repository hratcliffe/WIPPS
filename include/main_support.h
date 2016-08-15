//
//  main_support.h
//  
//
//  Created by Heather Ratcliffe on 25/02/2016.
//
//

#ifndef _main_support_h
#define _main_support_h


void get_deck_constants(std::string file_prefix);
int local_MPI_setup(int argc, char *argv[]);
setup_args process_command_line(int argc, char *argv[]);
std::vector<std::string> process_filelist(int argc, char *argv[]);
void share_consts();
void print_help(char code=0);
void divide_domain(std::vector<size_t>, int space[2], int per_proc, int block_num);

my_type get_ref_Bx(std::string file_prefix, int space_in[2], int time_0, bool is_acc=false);

#endif
