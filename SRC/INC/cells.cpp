#include                    "cells.h"
#include                    <cstdlib>
#include                    <cstring>
#include                    <cmath>
#include                    <iostream>
#include                    "inc.h"

void cells::init() {
  _neighbors = uptr_uset(new uset [POW_Nc[DIM]]);
  _num_of_neighbors = uptr_uint(new unsigned int [POW_Nc[DIM]]);
  _particles = uptr_uset(new uset [POW_Nc[DIM]]);
  _num_of_particles = uptr_uint(new unsigned int [POW_Nc[DIM]]);
  _cells_of_particles = uptr_uint(new unsigned int [N]);

  __loaded_once = false;

  FOR_PARTICLES(i) _cells_of_particles[i] = 0;
  FOR_CELLS(c) _num_of_particles[c] = 0;
  reload_neighbors();
  reload_particles();
}

/*
  This function builds the list of neighboring cells in which the periodic unit cell is divided.
  NOTE: if the shear (GAMMA) is 0 -> fine
        when GAMMA > 0, in order to account for all possible interactions at max distance = 2*R I
        have to add 3^(DIM-2) extra cells...

       #
    +--+--+--+
    |  #  |  |
  ==+==+==+==+== -> borders of the unit cell; neighbors span also in other unit cells (with periodic boundary conditions)
    |  # x|  | x and neighbors
    +--+--+--+
    |  #  |  |
    +--+--+--+
       #

    <->  #               <-> === GAMMA, translation
   +--+--+--+--+
   | y|  #  |  |   -> these cells are sheared
  ==+==+==+==+===
    |  # x|  |       x and neighbors; y is the extra neighbor, needed to cover all the interaction distance!
   -+--+--+--+--      (imagine the union of spheres of radius 2*R centered in every point in the cell x)
    |  #  |  |        (keep in mind that the side of each small cell is 1/N_c such that 1/N_c > 2*R)
   -+--+--+--+--
    |  #  |  |
*/
void cells::reload_neighbors() {
  unsigned int neighbor;

  if (GAMMA == 0) {
    FOR_CELLS(c) {
      _neighbors[c].clear();

      _num_of_neighbors[c] = POW_3[DIM];
      FOR_NUM_OF_NEIGHBORS(n,c) { // it refers to _num_of_neighbors just defined
        neighbor = 0;

        FOR_DIM(d) {
          unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
          int delta_digit = (int) (mod(floor(1.0*n/POW_3[d]), 3) - 1);
          int neighbor_digit = mod(int(cell_digit) + delta_digit, Nc);
          neighbor += neighbor_digit*POW_Nc[d];
        }

        _neighbors[c].insert(neighbor);
      }
    }

  } else {
    if (GAMMA < 0) std::cout<< "ERROR, GAMMA > 0!" << std::endl;////////////////////////////////////

    FOR_CELLS(c) {
      _neighbors[c].clear();
      unsigned int boundary_check = (unsigned int) mod(floor(1.0*c/Nc), Nc); // this is the 1-coordinate -> boundaries are at 0, Nc-1
      unsigned int shift = (unsigned int) floor(GAMMA*Nc);

      if (boundary_check == 0) {
        _num_of_neighbors[c] = POW_3[DIM];
        FOR_NUM_OF_NEIGHBORS(n,c) { // it refers to _num_of_neighbors just defined
          neighbor = 0;

          FOR_DIM(d) {
            unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
            int delta_digit = (int) (mod(floor(1.0*n/POW_3[d]), 3) - 1);
            int neighbor_digit = mod(int(cell_digit) + delta_digit, Nc);

            if ((neighbor_digit == int(Nc) - 1) && (d == 1)) { // in this case the new Nc-digit is 0, but I have to shift the 1-digit!
              neighbor_digit = mod(neighbor + shift, Nc);
              neighbor = neighbor_digit + POW_Nc[2] - Nc;
            } else {
              neighbor += neighbor_digit*POW_Nc[d];
            }
          }

          _neighbors[c].insert(neighbor);
        }

        _num_of_neighbors[c] += POW_3[DIM - 2];
        for (unsigned int counter = 0/*POW_3[DIM]*/; counter < /*POW_3[DIM] +*/ POW_3[DIM - 2]; counter++) { // extra neighbors
          neighbor = 0;

          int neighbor_digit = (int) mod(c + 2 + shift, Nc); // here c === its 1-digit
          neighbor += neighbor_digit + POW_Nc[2] - Nc;

          FOR_DIM(d) {
            if (d < 2) continue;
            unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
            int delta_digit = (int) (mod(floor(1.0*counter/POW_3[d]), 3) - 1);
            neighbor_digit = (int) mod(int(cell_digit) + delta_digit, Nc);

            neighbor += neighbor_digit*POW_Nc[d];
          }

          _neighbors[c].insert(neighbor);
        }

      } else if (boundary_check == Nc - 1) {
        _num_of_neighbors[c] = POW_3[DIM];
        FOR_NUM_OF_NEIGHBORS(n,c) { // in this case the new Nc-digit is 0, but I have to shift the 1-digit!
          neighbor = 0;

          FOR_DIM(d) {
            unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
            int delta_digit = (int) (mod(floor(1.0*n/POW_3[d]), 3) - 1);
            int neighbor_digit = (int) mod(int(cell_digit) + delta_digit, Nc);

            if ((neighbor_digit == 0) && (d == 1)) { // in this case the new Nc-digit is Nc-1, but I have to shift the 1-digit!
              neighbor_digit = (unsigned int) mod(int(neighbor) - shift, Nc);
              neighbor = neighbor_digit;
            } else {
              neighbor += neighbor_digit*POW_Nc[d];
            }
          }

          _neighbors[c].insert(neighbor);
        }

        _num_of_neighbors[c] += POW_3[DIM - 2];
        for (unsigned int counter = 0/*POW_3[DIM]*/; counter < /*POW_3[DIM] +*/ POW_3[DIM - 2]; counter++) { // extra neighbors
          neighbor = 0;

          int neighbor_digit = (int) mod(int(c) - 2 - shift, Nc); // here c === its 1-digit
          neighbor += neighbor_digit;

          FOR_DIM(d) {
            if (d < 2) continue;
            unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
            int delta_digit = (int) (mod(floor(1.0*counter/POW_3[d]), 3) - 1);
            neighbor_digit = (int) mod(int(cell_digit) + delta_digit, Nc);

            neighbor += neighbor_digit*POW_Nc[d];
          }

          _neighbors[c].insert(neighbor);
        }

      } else {
        _num_of_neighbors[c] = POW_3[DIM];
        FOR_NUM_OF_NEIGHBORS(n,c) {
          neighbor = 0;

          FOR_DIM(d) {
            unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
            int delta_digit = (int) (mod(floor(1.0*n/POW_3[d]), 3) - 1);
            int neighbor_digit = (int) mod(int(cell_digit) + delta_digit, Nc);

            neighbor += neighbor_digit*POW_Nc[d];
          }

          _neighbors[c].insert(neighbor);
        }
      }
    }
  }
}

/*cells& cells::operator= (const cells &copy) {
  if (this != &copy) {
    _neighbors = uptr_uset(new uset [POW_Nc[DIM]]);
    _num_of_neighbors = uptr_uint(new unsigned int [POW_Nc[DIM]]);
    _particles = uptr_uset(new uset [POW_Nc[DIM]]);
    _num_of_particles = uptr_uint(new unsigned int [POW_Nc[DIM]]);
    _cells_of_particles = uptr_uint(new unsigned int [N]);
    __loaded_once = copy.__loaded_once;
  }
  return *this;
}*/

void cells::reload_particles() {
  if (__loaded_once == false) {
    FOR_CELLS(c) _num_of_particles[c] = 0;

    FOR_PARTICLES(i) {
      unsigned int c_new = 0;
      FOR_DIM(d) c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];

      _cells_of_particles[i] = c_new;
      ++_num_of_particles[c_new];
      _particles[c_new].insert(i);
    }

    __loaded_once = true;

  } else {

    FOR_PARTICLES(i) {
      unsigned int c = cell_of(i), c_new = 0;
      FOR_DIM(d) c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];

      if (c != c_new) {
        _particles[c].erase(i);
        --_num_of_particles[c];

        _cells_of_particles[i] = c_new;
        ++_num_of_particles[c_new];
        _particles[c_new].insert(i);
      }
    }
  }

  return;
}

void cells::move_particle(unsigned int i, unsigned int c, unsigned int c_new) {
  _particles[c].erase(i);
  --_num_of_particles[c];

  _cells_of_particles[i] = c_new;
  ++_num_of_particles[c_new];
  _particles[c_new].insert(i);

  return;
}
