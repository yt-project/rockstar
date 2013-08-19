#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include "universal_constants.h"
#include "config_vars.h"
#include "universe_time.h"

#define STEPS 1024

static double t_0 = 0;       // Time now (in Hubble time units).
static double times[STEPS+1]={0};
static double H_CONV = HUBBLE_TIME_CONVERSION/0.7;
static double exact_t0_conv = 0;

void init_time_table(void) {
  double a = 1;
  int64_t i;
  H_CONV = HUBBLE_TIME_CONVERSION/h0;

  times[STEPS] = t_0;
  exact_t0_conv = exact_scale_to_time(1.0);
  for (i=1; i<=STEPS; i++) {
    a = 1.0-(double)i/((double)STEPS);
    times[STEPS-i] = exact_scale_to_time(a);
  }
}

double _exact_time_to_scale(double t) {
  double m = sinh(1.5*t*sqrt(1-Om));
  return pow(Om*m*m/(1.0-Om), 1.0/3.0);
}

double exact_scale_to_time(double scale) {
  double t = scale;
  double a = _exact_time_to_scale(t);
  double dt = scale/10.0;
  int64_t count = 0;
  while (fabs(a-scale)>1e-7 && (count < 10)) {
    count++;
    double a2 = _exact_time_to_scale(t+dt);
    double move = (scale-a)*(dt)/(a2-a);
    t += move;
    a = _exact_time_to_scale(t);
    if (move/10.0 < 0.5*dt) dt = move/10.0;
    else dt /= 2.0;
  }
  return t-exact_t0_conv;
}


// Linearly interpolate between calculated values.
double scale_to_time(double scale) {
  double s = scale;
  int64_t l = (int)(s*STEPS);
  double f = s*STEPS - l;
  if (scale > 1) return exact_scale_to_time(scale);
  if (scale < 0) return times[0];
  return (times[l]+f*(times[l+1]-times[l]));
}

double scale_to_years(double scale) {
  return (scale_to_time(scale)*H_CONV);
}


/*int main(void) {
  int64_t i;
  init_time_table();
  for (i=0; i<150; i++) {
    double a = _exact_time_to_scale(i/100.0);
    printf("%f %f %f %f\n", a, i/100.0-exact_t0_conv, exact_scale_to_time(a), scale_to_time(a));
  }
  }*/
