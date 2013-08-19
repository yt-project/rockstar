#ifndef HALO_H
#define HALO_H

#include <stdint.h>

struct halo {
  int64_t id;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    bullock_spin, kin_to_pot;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;
};

struct extra_halo_info {
  int64_t child, next_cochild, prev_cochild;
  int64_t sub_of, ph;
  float max_metric;
};

#endif /* HALO_H */
