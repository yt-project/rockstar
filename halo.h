#ifndef HALO_H
#define HALO_H

#include <stdint.h>

struct halo {
  int64_t id;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, mgrav, vmax, rvmax, rs, vrms, J[3], energy, spin;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;
};

struct extra_halo_info {
  int64_t child, next_cochild, prev_cochild;
  int64_t sub_of; //, cur_parent;
  float max_metric; //, prev_mass;
};

#endif /* HALO_H */
