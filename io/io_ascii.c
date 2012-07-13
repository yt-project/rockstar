#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "../check_syscalls.h"
#include "stringparse.h"
#include "../particle.h"

void gzip_file(char *filename) {
  char buffer[1024];
  snprintf(buffer, 1024, "gzip -f \"%s\"", filename);
  system(buffer);
}


void load_particles(char *filename, struct particle **p, int64_t *num_p) {
  FILE *input;
  char buffer[1024];
  int64_t n;
  struct particle d = {0};
  SHORT_PARSETYPE;
#define NUM_INPUTS 7
  enum short_parsetype stypes[NUM_INPUTS] = 
    { F, F, F, F, F, F, D64};
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(d.pos[0]), &(d.pos[1]), &(d.pos[2]), &(d.pos[3]), &(d.pos[4]), &(d.pos[5]), &(d.id)};

  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);

    if (n < NUM_INPUTS) continue;
    if (((*num_p)%1000)==0) {
      *p = check_realloc(*p, ((*num_p)+1000)*sizeof(struct particle),
			 "Adding new particles.");
    }
    (*p)[(*num_p)] = d;
    (*num_p)++;
  }
  fclose(input);
}
