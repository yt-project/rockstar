#ifndef _FUN_TIMES_H_
#define _FUN_TIMES_H_

#include <inttypes.h>

struct previous_halo {
  float pos[6];
  int64_t file_offset;
  int64_t p_offset;
  int64_t num_p;
  int64_t chunk;
  float m, r;
};

#define MAX_CORE_PARTICLES 10000
#define PREV_HALO_BUFFER_SIZE 100000

void load_previous_halos(int64_t snap, int64_t chunk, float *bounds);
float find_previous_mass(struct halo *h, struct particle *hp, int64_t *best_num_p);
void convert_and_sort_core_particles(struct halo *h, struct particle *hp);

#endif /* _FUN_TIMES_H_ */
