#include <stdio.h>
#include <math.h>
#include "universal_constants.h"

inline double c_to_f(double c) {
  double cp1 = 1.0+c;
  return (c*cp1 / (log1p(c)*cp1 - c));
}

double f_to_c(double f) {
  double c = f;
  double tc = c_to_f(c);
  double new_c;
  while (fabs((f-tc) / f) > 1e-7) {
    double tc2 = c_to_f(c+0.1);
    double slope = (tc2 - tc) / 0.1;
    new_c = c + (f-tc)/slope;
    if (new_c < 0) c/=2;
    else c = new_c;
    tc = c_to_f(c);
  }
  return c;
}

float calc_scale_radius(float mvir, float rvir, float vmax, float rvmax, float scale)
{
  float f, c, vm2;
  if (!mvir || !rvir || !vmax) return (rvmax / RMAX_TO_RS);
  vm2 = vmax/VMAX_CONST;
  vm2 = vm2*vm2 * scale;
  f = (rvir/1.0e3)*vm2 / (mvir*RS_CONSTANT);
  if (f < 4.625) return (rvmax / RMAX_TO_RS);
  c = f_to_c(f);
  if (c <= 0) return (rvmax / RMAX_TO_RS);
  return (rvir/c);
}
