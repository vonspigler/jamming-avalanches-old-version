#ifndef                     HEADERS_FUNCTS
#define                     HEADERS_FUNCTS

#include                    <string>

#define FOR_CNT(cnt)\
  for (auto __local_counter__ = 0; __local_counter__ < cnt; ++__local_counter__)

/*##################################################################################################
##                                                                                                ##
##   - minimize()                                                                                 ##
##   - find_jamming()                                                                             ##
##   - delta_shear() shears the system with some given parameter                                  ##
##   - print_conf() print current configuration in a separate file (const path...)                ##
##   - check_args() controls **argv and sets some values (N,DIM)                                  ##
##   - init_all(), init_phi_pos() perform some initilizations at different times (phi>0?)         ##
##   - print_minimization_time() cout's the time for a minimization run (also minimize()'s)       ##
##   - adjust_parameters() adjusts R,Nc according to PHI                                          ##
##                                                                                                ##
##################################################################################################*/

bool minimize(); // returns true if it took ITER_MAX steps to finish the minimization
void find_jamming();
void print_conf(std::string);
int check_args(int, char**);
void init_all();

#endif
