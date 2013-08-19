#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <stdint.h>
#define POTENTIAL_DONT_CALCULATE_FLAG 1

struct potential {
  float pos[6], r2;
  double pe;
  float ke;
#ifdef POTENTIAL_COMPARISON
  float pe2, ke2, pe3, ke3;
  float v,r;
#endif /* POTENTIAL_COMPARISON */
#ifdef CALC_POTENTIALS
  int64_t id;
#endif /* CALC_POTENTIALS */
  int32_t flags;
};

void compute_kinetic_energy(struct potential *po, int64_t num_po, float *vel_cen, float *pos_cen);
void compute_potential(struct potential *po, int64_t num_po);

#endif /* _POTENTIAL_H_ */
