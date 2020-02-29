#ifndef                     HEADERS_INC
#define                     HEADERS_INC

#include                    <memory> // c++11
#include                    <string>
#include                    <random>
#include                    <vector>
#include                    "vectors.h"
#include                    "cells.h"

#define FOR_CELLS(c)        for (unsigned int c = 0; c < POW_Nc[DIM]; ++c)
#define FOR_NEIGHBORS(x,c)  for (const auto x : CELL_LIST.neighbors(c))
#define FOR_DIM(x)          for (unsigned int x = 0; x < DIM; ++x)
#define FOR_PARTICLES(x)    for (unsigned int x = 0; x < N; ++x)
#define FOR_PARTICLES_IN_CELL(x,c)\
                            for (const auto x : CELL_LIST.particles(c)
#define FOR_PARTICLES_NEAR_CELL(x,c)\
                            for (const auto n : CELL_LIST.neighbors(c)) for (const auto x : CELL_LIST.particles(n))
#define FOR_PARTICLES_NEAR_PARTICLE(x,y)\
                            for (const auto n : CELL_LIST.neighbors(CELL_LIST.cell_of(y))) for (const auto x : CELL_LIST.particles(n))

/*##################################################################################################
##                                                                                                ##
##   - extern variables (R,N,DIM,PHI,Nc,precisions,GAMMA,ENERGY,XS,YS,AS,CELL_LIST)               ##
##   - extern MD/FIRE consts                                                                      ##
##   - extern POWers (3,Nc)                                                                       ##
##   - MPI references to proc. num.                                                               ##
##   - random wrappers                                                                            ##
##   - mod function                                                                               ##
##   - energy functions, derivatives                                                              ##
##   - contacts, overlap, distances functions                                                     ##
##                                                                                                ##
##################################################################################################*/

extern unsigned int N;
extern unsigned int Nc;
extern unsigned int DIM;
extern double R;
extern const double PREC;
extern const double PREC_PHI;
extern const double PREC_Z;
extern const unsigned int ITER_MAX;
extern double GAMMA;
extern double PHI;

extern std::string OUTPUT_FOLDER;

extern std::unique_ptr<unsigned int []> POW_3;
extern std::unique_ptr<unsigned int []> POW_Nc;

extern const double MD_alpha_start;
extern const double MD_delta_t_max;
extern const double MD_f_alpha;
extern const double MD_f_inc;
extern const double MD_f_dec;
extern const unsigned int MD_N_min;
extern double MD_alpha;
extern double MD_delta_t;
extern unsigned int MD_counter;
extern double MD_V_norm;
extern double MD_A_norm;

extern double __dim_param;

extern int __world_rank;
extern int __world_size;

extern std::default_random_engine __RAND_GENERATOR;
extern std::uniform_int_distribution<unsigned int> __RAND_INT;
extern std::uniform_real_distribution<double> __RAND_DOUBLE;
extern std::normal_distribution<double> __RAND_GAUSS;

inline unsigned int rand_int(unsigned int n)
  { return (unsigned int)(1.0*n*__RAND_INT(__RAND_GENERATOR)/__RAND_INT.max()); }
inline double rand_double(double n)
  { return n*__RAND_DOUBLE(__RAND_GENERATOR); }
inline double rand_gauss(double m, double s)
  { return __RAND_GAUSS(__RAND_GENERATOR)*s + m; }
inline double abs_v(double x) { return (x < 0 ? -x : x); }

std::string get_seed();
void set_seed(std::string);

inline bool draw_prob(double p) { return (1.0*rand()/RAND_MAX < p); }
inline double mod(double x, unsigned int n) { return x - n*floor(x/n); }

//double H(); // total energy
//double H(unsigned int); // energy felt by one particle
//double H(unsigned int, unsigned int); // couple interaction between two particles
//double H_prime(unsigned int, unsigned int); // force component on particle
void compute_forces(vectors<double> &);

double dist_1d_pbcs(unsigned int, unsigned int, const unsigned int); // corrected pbcs distance along dimension d
double dist_1d_pbcs(unsigned int, unsigned int, const unsigned int, const vectors<double> &); // same as before, but between XS and given vector
double sq_dist_pbcs(unsigned int, unsigned int); // periodic boundary conditions
double sq_dist_pbcs(unsigned int, unsigned int, const vectors<double> &); // same as before, but between XS and given vector

//double force_between_particles(unsigned int, unsigned int, unsigned int);
double contacts_avg_pp();
void compute_observables(double &, double &, double &, double &, const vectors<double> &); // energy, contacts, overlap, stress

//void perm(std::vector<unsigned int> &, unsigned int);
void MDFIRE_step();

extern vectors<double> XS;
extern vectors<double> VS;
extern vectors<double> AS;
extern cells CELL_LIST;

#endif  
