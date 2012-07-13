#ifndef _BOUNDS_H_
#define _BOUNDS_H_

int _check_bounds(float *pos_i, float *pos_f, float *bounds);
int _check_bounds_raw(float *pos_i, float *bounds);
int bounds_overlap(float *b1, float *b2, float *b3, double overlap);

#endif /* _BOUNDS_H_ */
