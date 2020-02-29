#include                    "other.h"
#include                    <cmath>
//#include                    <chrono> // c++11
#include                    <iostream>
#include                    <fstream>
#include                    <cstdlib>
#include                    <memory> // c++11
#include                    <string>
#include                    "inc.h"

void find_jamming() {
  double PHI_0 = 0.0, PHI_1 = 1.0, z_num = 0;

  do {
    PHI = (PHI_0 + PHI_1)/2;
    init_all(); // minimizes too!

    z_num = contacts_avg_pp();
//std::cout<<R<<" "<<PHI<<" "<<z_num<<std::endl;
//print_conf();

    if (z_num > 2*DIM) PHI_1 = PHI;
    else PHI_0 = PHI;
std::cout << fabs(z_num - 2*DIM) << std::endl;
  } while ((fabs(PHI_1 - PHI_0) > PREC_PHI) || (fabs(z_num - 2*DIM) > PREC_Z));

  PHI += PREC_Z/10; 
  init_all();

// Now I'm near the jamming point. z-z_c ~ sqrt[(phi-phi_c)^+]; find phi_c fitting (z-z_c)^2 to phi-phi_c
// PROBLEM: the sqrt ends in (~0.639,5.82)!!! z_c != 6.... with 1000 particles!?
// -> RATTLERS! I should remove them with some kind of leaf removal algorithm!
////////////////////////////////////////////////////////////////////////////////////////////////////
//  if (false) {                                                                                  //
//    double old_phi = PHI;                                                                       //
                                                                                                  //
//    for (double pckg_f = old_phi - 0.01; pckg_f <= old_phi + 0.1; pckg_f += 0.005) {            //
//      PHI = pckg_f;                                                                             //
//      R = pow(PHI/N*param, 1.0/DIM);                                                            //
//      Nc = ceil(1./2/R - 1);                                                                    //
//      if (Nc < 1) {                                                                             //
//        std::cout << "Nc >= 1 && R <= 1/2!" << std::endl;                                       //
//        return;                                                                                 //
//      }                                                                                         //
                                                                                                  //
//      CELL_LIST.reload_particles();                                                             //
//      CELL_LIST.reload_neighbors();                                                             //
//      compute_forces(AS)                                                                        //
                                                                                                  //
//      minimize();                                                                               //
//      std::cout << PHI <<" "<< contacts_avg_pp() << std::endl;                                  //
//    }                                                                                           //
                                                                                                  //
//    PHI = old_phi;                                                                              //
//  }                                                                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
  return;
}

void init_all() {
  R = pow(PHI/N*__dim_param, 1.0/DIM);
  Nc = ceil(1./2/R - 1);
  if (Nc < 1) std::cout << "Nc >= 1 && R <= 1/2!" << std::endl;

  POW_3 = std::unique_ptr<unsigned int []>(new unsigned int [DIM + 1]);
  POW_Nc = std::unique_ptr<unsigned int []>(new unsigned int [DIM + 1]);
  POW_3[0] = 1;
  POW_Nc[0] = 1;

  for (unsigned int d = 1; d < DIM + 2; d++) {
    POW_3[d] = 3*POW_3[d - 1];
    POW_Nc[d] = Nc*POW_Nc[d - 1];
  }

  CELL_LIST = cells(); // necessary?
  CELL_LIST.init();
  compute_forces(AS);
  minimize();

 return;
}

int check_args(int argc, char **argv) { // also computes N,DIM
  if (argc != 9) {
    std::cout << "Use: ./MAIN.bin N DIM Phi SAMPLES_per_processor Delta_GAMMA GAMMA_MAX SEED OUTPUT_FOLDER/" << std::endl;
    return -1;
  }

  N = atoi(argv[1]);
  if (N < 1) {
    std::cout << "N >= 1!" << std::endl;
    return -1;
  }

  DIM = atoi(argv[2]);
  if (DIM <= 1) {
    std::cout << "DIM > 1!" << std::endl;
    return -2;
  }

  OUTPUT_FOLDER = std::string(argv[8]);

  return 0;
}

bool minimize() {
  MD_alpha = MD_alpha_start;
  MD_delta_t = MD_delta_t_max/10;
  MD_counter = 0;
  unsigned int counter = 0;

  do {
    MDFIRE_step();
  } while ((++counter < ITER_MAX*N) && (MD_A_norm > PREC));

  if (counter == ITER_MAX*N) return false;
  return true;
}

void print_conf(std::string name) {
  std::ofstream OUT_conf;
  std::string n(OUTPUT_FOLDER);
  n += "conf";
  n += name;
  OUT_conf.open(n.c_str(), std::ios::out);

  FOR_PARTICLES(i) {
    OUT_conf << R << " ";
    FOR_DIM(d) OUT_conf << XS(i,d) << " ";
    OUT_conf << std::endl;
  }

  OUT_conf.close();

  return;
}

/*
auto __start_time = std::chrono::system_clock::now();
auto __stop_time = std::chrono::system_clock::now();
std::cout << " -- " << ((std::chrono::duration<double>)(__stop_time - __start_time)).count() << "s\n";
*/
