#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <stdint.h>
#define POTENTIAL_DONT_CALCULATE_FLAG 1

struct potential {
  float pos[6], r2;
  //union {
    double pe;
  //  int64_t haloid;
  //};
  //union {
    float ke;
  //  float density;
  //};
#ifdef POTENTIAL_COMPARISON
  double pe2;
  float v,r;
#endif /* POTENTIAL_COMPARISON */
  int32_t flags;
};

void compute_kinetic_energy(struct potential *po, int64_t num_po, float *vel_cen, float *pos_cen);
void compute_potential(struct potential *po, int64_t num_po);

#endif /* _POTENTIAL_H_ */
