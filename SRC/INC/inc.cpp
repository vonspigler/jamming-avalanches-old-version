#include                    "inc.h"
#include                    <cstring>
#include                    <cstdlib>
#include                    <cmath>
#include                    <string>
#include                    <sstream>
#include                    <vector>
#include                    <iostream>
//#include                    <chrono> // c++11

#define                     WITH_PROB(p)    if (draw_prob(p))

/*-------  RANDOM NUMBERS  -----------------------------------------------------------------------*/

std::string get_seed() {
  std::stringbuf s;
  std::ostream __out(&s);
  __out << __RAND_GENERATOR;
  return s.str();
}

void set_seed(std::string gen) {
  std::stringbuf s;
  s.str(gen);
  std::istream __int(&s);
  __int >> __RAND_GENERATOR;

  __RAND_INT.reset();
  __RAND_DOUBLE.reset();
  __RAND_GAUSS.reset();

  return;
}

double sq_dist_pbcs(unsigned int i, unsigned int j) {
  double dist, r = 0, del;

  short n;
  del = XS(j,1) - XS(i,1);

  if (del > 0.5) n = 1;
  else if (del < -0.5) n = -1;
  else n = 0;

  del = XS(j,0) - XS(i,0) - GAMMA*n;

  r += del;
  if (del > 0.5) --r;
  else if (del < -0.5) ++r;

  dist = r*r;

  for (unsigned int d = 1; d < DIM; ++d) {  // only d > 0
    r = 0;

    del = XS(j,d) - XS(i,d);

    r += del;
    if (del > 0.5) --r;
    else if (del < -0.5) ++r;

    dist += r*r;
  }

  return dist;
}

double sq_dist_pbcs(unsigned int i, unsigned int j, const vectors<double> &YS) {
  double dist, r = 0, del;

  short n;
  del = XS(j,1) - YS(i,1);

  if (del > 0.5) n = 1;
  else if (del < -0.5) n = -1;
  else n = 0;

  del = XS(j,0) - YS(i,0) - GAMMA*n;

  r += del;
  if (del > 0.5) --r;
  else if (del < -0.5) ++r;

  dist = r*r;

  for (unsigned int d = 1; d < DIM; ++d) {  // only d > 0
    r = 0;

    del = XS(j,d) - YS(i,d);

    r += del;
    if (del > 0.5) --r;
    else if (del < -0.5) ++r;

    dist += r*r;
  }

  return dist;
}

double dist_1d_pbcs(unsigned int j, unsigned int i, const unsigned int d) {
  double del;

  if (d != 0) {
    del = XS(j,d) - XS(i,d);
    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;


  } else if (d == 0) {
    short n;
    del = XS(j,1) - XS(i,1);

    // I AM ACTUALLY ASSUMING THAT |n| <= 1 --> ok since the two particles will be belonging to neighboring cells
    if (del > 0.5) n = 1;
    else if (del < -0.5) n = -1;
    else n = 0;

    del = XS(j,0) - XS(i,0) - GAMMA*n;

    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;
  }

  return -1;
}

double dist_1d_pbcs(unsigned int j, unsigned int i, const unsigned int d, const vectors<double> &YS) {
  double del;

  if (d != 0) {
    del = XS(j,d) - YS(i,d);
    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;

  } else if (d == 0) {
    short n;
    del = XS(j,1) - YS(i,1);

    if (del > 0.5) n = 1;
    else if (del < -0.5) n = -1;
    else n = 0;

    del = XS(j,0) - YS(i,0) - GAMMA*n;

    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;
  }

  return -1;
}

/*-------  ENERGY FUNCTION  -------------------*/

// I am actually not using those
/*
double H() {
  double r = 0, dist;

  FOR_PARTICLES(i) FOR_PARTICLES_NEAR_PARTICLE(j,i) {
    if (i != j) {
      dist = sqrt(sq_dist_pbcs(i, j));
      if (dist < 2*R) r += (1 - dist/2/R)*(1 - dist/2/R);
    }
  }

  return r/2;
}

double H(unsigned int i, unsigned int j) {
  double r = 0, dist;
  dist = sqrt(sq_dist_pbcs(i, j));
  if (dist < 2*R) r = (1 - dist/2/R)*(1 - dist/2/R);
  return r;
}

double H(unsigned int i) {
  double r = 0, dist;
  unsigned int c = CELL_LIST.cell_of(i);

  FOR_PARTICLES_NEAR_CELL(j,c) {
    if (i != j) {
      dist = sqrt(sq_dist_pbcs(i, j));
      if (dist < 2*R) r += (1 - dist/2/R)*(1 - dist/2/R);
    }
  }

  return r;
}


double H_prime(unsigned int i, unsigned int d) {
  double r = 0, dist;

  FOR_PARTICLES_NEAR_PARTICLE(j,i) {
    if (i != j) {
      dist = sqrt(sq_dist_pbcs(i, j));
      if (dist < 2*R) r -= (1 - dist/2/R)*dist_1d_pbcs(i, j, d)/R/dist;
    }
  }

  return r;
}*/

void compute_forces(vectors<double> &a) {
  double r, dist;

  FOR_PARTICLES(i) FOR_DIM(d) {
    r = 0;

    FOR_PARTICLES_NEAR_PARTICLE(j,i) {
      if (i != j) {
        dist = sqrt(sq_dist_pbcs(i, j));
        if (dist < 2*R) r += (1 - dist/2/R)*dist_1d_pbcs(i, j, d)/R/dist;  // force = -h_prime!
      }
    }

    a(i,d) = r;
  }

  return;
}

// this could be useful if opening contacts rather than shearing
/*
double force_between_particles(unsigned int i, unsigned int j, unsigned int d) {
  double r = -1, dist;
  if (i != j) {
    dist = sqrt(sq_dist_pbcs(i, j));
    if (dist < 2*R) r = -(1 - dist/2/R)*dist_1d_pbcs(i, j, d)/R/dist;
  }

  return r;
}
*/

double contacts_avg_pp() {
  double z = 0;
  double dist_sq;

  FOR_PARTICLES(i) FOR_PARTICLES_NEAR_PARTICLE(j,i) {
    if (i != j) {
      dist_sq = sq_dist_pbcs(i, j);
      if (dist_sq < 4*R*R) z += 1;
    }
  }

  return z/N;
}

void compute_observables(double &e, double &z, double &o, double &syx, const vectors<double> &YS) {
  double t_e = 0, t_z = 0, t_o = 0, t_s_yx = 0/*, t_s_xy = 0*/, dist_sq, dist;

  FOR_PARTICLES(i) FOR_PARTICLES_NEAR_PARTICLE(j,i) {
    if (i != j) {
      dist_sq = sq_dist_pbcs(i, j);

      if (dist_sq < 4.0*R*R) {
        dist =sqrt(dist_sq);

        t_e += (1 - dist/2/R)*(1 - dist/2/R);
        t_z += 1;
        t_s_yx -= (1 - dist/2/R)*dist_1d_pbcs(i, j, 0)*dist_1d_pbcs(i, j, 1)/R/dist;
      }
    }

    dist_sq = sq_dist_pbcs(i, j, YS);
    if (dist_sq < 4.0*R*R) t_o += exp(-dist_sq/R/R);
  }

  e = t_e/2;
  z = t_z/N;
  o = t_o/N;
  syx = t_s_yx/2/N;
  return;
}

/*void perm(std::vector<unsigned int> &s, unsigned int max) {
  unsigned int counter = max, k, tmp;

  do {
    k = rand_int(counter - 2);
    tmp = s[k];
    s[k] = s[counter - 1];
    s[counter - 1] = tmp;
    --counter;
  } while(counter > 1);

  return;
}*/

// forces should be computed before calling this function
void MDFIRE_step() {
  int n_disp;

  // Velocity Verlet
  FOR_PARTICLES(i) {
    // The y-component has to be updated before the x-component!!! ???
    FOR_DIM(d) {
      XS(i,d) += VS(i,d)*MD_delta_t + AS(i,d)*MD_delta_t*MD_delta_t/2;
      VS(i,d) += AS(i,d)*MD_delta_t/2; // updating V's and X's
    }

    n_disp = floor(XS(i,1));
    if (n_disp != 0) XS(i,0) -= GAMMA*n_disp; // no VS, AQS

    FOR_DIM(d) if ((XS(i,d) >= 1) || (XS(i,d) < 0)) XS(i,d) -= floor(XS(i,d));

    unsigned int c_new = 0, tmp_c = CELL_LIST.cell_of(i);
    FOR_DIM(d)  c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];
    if (tmp_c != c_new) CELL_LIST.move_particle(i, tmp_c, c_new);
  }

  compute_forces(AS);
  MD_A_norm = 0;
  MD_V_norm = 0;

  FOR_PARTICLES(i) FOR_DIM(d) {
    VS(i,d) += AS(i,d)*MD_delta_t/2;

    MD_V_norm += VS(i,d)*VS(i,d);
    MD_A_norm += AS(i,d)*AS(i,d);
  }

  MD_V_norm = sqrt(MD_V_norm);
  MD_A_norm = sqrt(MD_A_norm);

  double A_norm, P = 0;

  // FIRE
  FOR_PARTICLES(i) FOR_DIM(d) {
    A_norm = (MD_A_norm >= 1e-8 ? AS(i,d)/MD_A_norm : 0);

    VS(i,d) = (1 - MD_alpha)*VS(i,d) + MD_alpha*MD_V_norm*A_norm;
    P += VS(i,d)*AS(i,d);
  }

  if (P <= 0) {
    FOR_PARTICLES(i) FOR_DIM(d) VS(i,d) = 0;

    MD_alpha = MD_alpha_start;
    MD_counter = 0;
    MD_delta_t *= MD_f_dec;

  } else {
    if (MD_counter > MD_N_min) {
      MD_delta_t *= MD_f_inc;
      if (MD_delta_t > MD_delta_t_max) MD_delta_t = MD_delta_t_max;

      MD_alpha *= MD_f_alpha;
    }

    MD_counter++;
  }

/////////////////////////////////////////////////////////////////////////////////////
/*double energy = 0, dist, cz = 0;
FOR_PARTICLES(i) FOR_PARTICLES_NEAR_PARTICLE(j,i) {
  if (i != j) {
    dist = sqrt(sq_dist_pbcs(i, j));
    if (dist < 2*R) {
      energy += (1 - dist/2/R)*(1 - dist/2/R);
      ++cz;
    }
  }
}
std::cout << MD_A_norm << " " << energy/2 << " " << cz/N << "\n";*/
///////////////////////////////////////////////////////////////////////////////////

  return;
}
