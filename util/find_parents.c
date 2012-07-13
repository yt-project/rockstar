#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../check_syscalls.h"
#include "../io/stringparse.h"
#include <math.h>

#define EXTRA_HALO_INFO int64_t descid, np;
#include "read_tree.h"

double BOX_SIZE=250;

struct halo_list all_halos = {0};
struct halo_tree halo_tree = {0};

#define GROUP_LIST all_halos.halos
#define RADIUS rvir
#define FAST3TREE_TYPE struct halo
#include "../fast3tree.c"
#define parent pid
#include "../parents.c"
#undef parent

void read_hlist(char *filename) {
  int64_t n;
  FILE *input;
  struct halo h = {0};
  char buffer[1024];

  SHORT_PARSETYPE;
  #define NUM_INPUTS 14
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F, F, F,    //  #id desc_id mvir vmax vrms
      F, F, D64, F,       //  Rvir Rs Np x y z vx vy vz 
      F, F, F, F, F,    
    };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(h.id),
                            &(h.descid),			    
			    &(h.mvir), &(h.vmax), &(h.vrms), &(h.rvir), &(h.rs), 
			    &(h.np), 
			    &(h.pos[0]), &(h.pos[1]), &(h.pos[2]), 
			    &(h.vel[0]), &(h.vel[1]), &(h.vel[2])};
  

  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    if (n<NUM_INPUTS) continue;
    if (!(all_halos.num_halos%3000))
      all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*(all_halos.num_halos+3000), "Allocating Halos.");
   
    all_halos.halos[all_halos.num_halos] = h;
    all_halos.num_halos++;
  }
  fclose(input);
  
  all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*all_halos.num_halos, "Allocating Halos.");

  find_parents(all_halos.num_halos);

  printf("#ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ PID\n");
  for (n=0; n<all_halos.num_halos; n++) {
    struct halo *th = all_halos.halos + n;
    printf("%"PRId64" %"PRId64" %.3e %.2f %.2f %.3f %.3f %"PRId64" %.5f %.5f %.5f %.2f %.2f %.2f %"PRId64"\n",
	  th->id, th->descid, th->mvir, th->vmax, th->vrms, th->rvir, th->rs, th->np,
	  th->pos[0], th->pos[1], th->pos[2], th->vel[0], th->vel[1],
	   th->vel[2], th->pid);
  }
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    printf("Usage: %s hlist box_size\n", argv[0]);
    exit(1);
  }
  if (argc > 2) BOX_SIZE = atof(argv[2]);
  read_hlist(argv[1]);
  return 0;
}


