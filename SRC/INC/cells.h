#ifndef                     HEADERS_CELLS
#define                     HEADERS_CELLS

#include                    <unordered_set> // c++11
#include                    <memory>

#define FOR_NUM_OF_NEIGHBORS(x,c)\
  for (unsigned int x = 0; x < CELL_LIST.num_of_neighbors(c); x++)

/*##################################################################################################
##                                                                                                ##
##   - class CELLS                                                                                ##
##                                                                                                ##
##   It organizes the subcells of the simulation, loading the particles etc...                    ##
##                                                                                                ##
##################################################################################################*/

using uset = std::unordered_set<unsigned int>;
using uptr_uset = std::unique_ptr<uset []>;
using uptr_uint = std::unique_ptr<unsigned int []>;

class cells {
private:
  uptr_uset _neighbors;
  uptr_uint _num_of_neighbors;
  uptr_uset _particles;
  uptr_uint _num_of_particles;
  uptr_uint _cells_of_particles;
  bool __loaded_once;

public:
  void init();
  void reload_particles();
  void move_particle(unsigned int, unsigned int, unsigned int);
  void reload_neighbors();

  inline const uset &neighbors(unsigned int c) const { return _neighbors[c]; };
  inline const uset &particles(unsigned int c) const { return _particles[c]; };
  inline unsigned int cell_of(unsigned int i) const { return _cells_of_particles[i]; };
  inline unsigned int num_of_particles(unsigned int c) const { return _num_of_particles[c]; };
  inline unsigned int num_of_neighbors(unsigned int c) const { return _num_of_neighbors[c]; };
};

#endif
