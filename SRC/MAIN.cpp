#include                    <cmath>
#include                    <cstdlib>
#include                    <iostream>
#include                    <fstream>
#include                    <sstream>
#include                    <chrono> // c++11
#include                    <string>
#include                    <memory> // c++11
#include                    <mpi.h>
#include                    "INC/vectors.h"
#include                    "INC/cells.h"
#include                    "INC/inc.h"
#include                    "INC/other.h"

#define                     MASTER          0

/*####  PROGRAM NOTES  #############################################################################
##                                                                                                ##
##  - c++11 libraries: memory; unordered_set; chrono                                              ##
##                                                                                                ##
##  TODO:                                                                                         ##
##    constant pressure, 10^-5 (increasing volume, constant (or min?) H-pV                        ##
##    bigger cells or rather neighbor lists                                                       ##
##    statistics of particles and modes in avalanches                                             ##
##    play with sigma_f (PREC) and with sigma_e (induced) to justify noise                        ##
##    yielding                                                                                    ##
##    play with dg (scalings) and, perhaps, g_max                                                 ##
##    play with mass/dt/energy scale in H in order to minimize quickly                            ##
##    perhaps only save avalanches (eg displacement or energy threshold?)                         ##
##                                                                                                ##
######  PHYSICS NOTES  #############################################################################
##                                                                                                ##
##  OUTPUT:                                                                                       ##
##      gamma (accumulated), energy (=>total), z, overlap (t-1,t: bc non additive), stress xy     ##
##          energy and stress are intensive (change this, I want them to scale with N)            ##
##          I should also compute statistics of particles involved in avalances?                  ##
##          And eigenmodes/values                                                                 ##
##                                                                                                ##
##################################################################################################*/

using namespace             std;

/*######################################################################  GLOBAL VARIABLES  ######*/

unsigned int                N;
unsigned int                Nc;
unsigned int                DIM;
double                      R;
double                      GAMMA;
double                      PHI;

uptr_uint                   POW_3;
uptr_uint                   POW_Nc;

const unsigned int          ITER_MAX        = 10000;  // maximum number of cycles per particle for minimize()
const double                PREC            = 1e-4; // minimize() converges when |VS|, |AS| <= PREC
const double                PREC_PHI        = 1e-4;
const double                PREC_Z          = 1e-4;

default_random_engine                       __RAND_GENERATOR;
uniform_int_distribution<unsigned int>      __RAND_INT;
uniform_real_distribution<double>           __RAND_DOUBLE(0, 1);
normal_distribution<double>                 __RAND_GAUSS(0, 1);

std::string                 OUTPUT_FOLDER;

const double                MD_alpha_start  = 0.1;
const double                MD_delta_t_max  = 0.01; // check!
const double                MD_f_alpha      = 0.99;
const double                MD_f_inc        = 1.1;
const double                MD_f_dec        = 0.5;
const unsigned int          MD_N_min        = 5;
double                      MD_alpha;
double                      MD_delta_t;
unsigned int                MD_counter;
double                      MD_V_norm;
double                      MD_A_norm;

double                      __dim_param;

vectors<double>             XS;
vectors<double>             VS;
vectors<double>             AS;
cells                       CELL_LIST;

int                         __world_rank;
int                         __world_size;

/*##################################################################################  MAIN  ######*/

int main(int argc, char **argv) {
  if (check_args(argc, argv) < 0) return -1;

  MPI_Init(NULL, NULL); // different processes will run different samples
  MPI_Comm_size(MPI_COMM_WORLD, &__world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &__world_rank);

  unsigned int SAMPLES = atoi(argv[4]);
  double Delta_GAMMA = atof(argv[5]);
  double GAMMA_MAX = atof(argv[6]);

  int SEED = atoi(argv[7]);
  if (SEED > 0) {
    set_seed(to_string(__world_rank + SEED));
  } else {
    SEED = chrono::system_clock::now().time_since_epoch().count();
    set_seed(to_string(__world_rank + SEED));
  }

  vectors<double> XS_0(N, DIM), XS_prev(N,DIM);
  XS = vectors<double>(N, DIM);
  VS = vectors<double>(N, DIM);
  AS = vectors<double>(N, DIM);

  __dim_param = tgamma(DIM/2.0 + 1)/pow(M_PI, DIM/2.0);

  PHI = atof(argv[3]);
  double PHI_start = PHI;

  /*######################################################################  SAMPLE AVERAGE  ######*/

  for (unsigned int sample = 0; sample < SAMPLES; sample++) {
    FOR_PARTICLES(i) FOR_DIM(d) {
      XS(i,d) = rand_double(1.0);
      VS(i,d) = 0;
    }

    PHI = PHI_start;
    GAMMA = 0;
    if (PHI <= 0) find_jamming(); // they do also
    else init_all();              // minimize!

    FOR_PARTICLES(i) FOR_DIM(d) XS_0(i,d) = XS(i,d);

    double energy, contacts_num, overlap, stress_yx;
    compute_observables(energy, contacts_num, overlap, stress_yx, XS_0);

    auto n = OUTPUT_FOLDER + to_string(__world_rank*SAMPLES + sample);
    ofstream OUT;
    OUT.open(n.c_str(), ios::out);
//    print_conf(to_string(__world_rank*SAMPLES + sample) + "_init");

    OUT << "# PHI=" << PHI << ", R=" << R << ", Nc=" << Nc << ", SEED=" << SEED <<  endl;
    OUT << 0 << " " << energy << " " << contacts_num << " " << overlap << " " << stress_yx << endl;

    /*########################################################################  SHEAR LOOP  ######*/

//GAMMA = Delta_GAMMA; //////////////////////////////////////////// for small GAMMA_MAX (<1/Nc) I need
//CELL_LIST.reload_neighbors(); /////////////////////////////////// to reload_neighbors() only once!

    for (auto GAMMA_loop = Delta_GAMMA; GAMMA_loop <= GAMMA_MAX; GAMMA_loop += Delta_GAMMA) {
      GAMMA = GAMMA_loop;
      CELL_LIST.reload_neighbors();
      FOR_PARTICLES(i) {
        FOR_DIM(d) XS_prev(i,d) = XS(i,d);  // Jumps in overlap between t and t-1! Not t and 0.

        // Affine transformation
        XS(i,0) += Delta_GAMMA*XS(i,1); // no VS, AQS

        int n_disp = floor(XS(i,1));
        if (n_disp != 0) XS(i,0) -= GAMMA*n_disp; // no VS, AQS
        if ((XS(i,0) >= 1) || (XS(i,0) < 0)) XS(i,0) -= floor(XS(i,0));
      }

      CELL_LIST.reload_particles(); // I have updated the coordinates

      compute_forces(AS); // important
      if (!minimize()) cout << "# GAMMA=" << GAMMA_loop << ", FILE=" << n.c_str() << ", NOT CONVERGED" << endl;
      compute_observables(energy, contacts_num, overlap, stress_yx, XS_prev);
      OUT << GAMMA_loop << " " << energy << " " << contacts_num << " " << overlap << " " << stress_yx << endl;
    }


//    print_conf(to_string(__world_rank*SAMPLES + sample) + "_shear");
    OUT.close();
  }

  MPI_Finalize();
  return 0;
}
