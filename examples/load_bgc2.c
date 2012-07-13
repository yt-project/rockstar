/* Example code for loading BGC2 output files. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include "../check_syscalls.h"
#include "../io/bgc2.h"
#include "../io/io_util.h"

GROUP_DATA_RMPVMAX *grps = NULL;
PARTICLE_DATA_PV *parts = NULL;
int64_t num_g=0, num_p=0;

void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PV **pdata, int64_t *num_parts);


int main(int argc, char **argv) {
  int64_t i;
  struct bgc2_header hdr;

  if (argc < 2) {
    printf("Usage: %s file1.bgc2 ...\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc; i++) {
    num_g = num_p = 0;
    load_bgc2(argv[i], &hdr, &grps, &num_g, &parts, &num_p);
    //Do stuff with particles here...
    //See io/bgc2.h for a description of the group and particle structure format
  }
  return 0;
}



void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PV **pdata, int64_t *num_parts)
{
  FILE *input;
  int64_t new_group_size, new_part_size;
  int64_t i, p_start;

  assert(sizeof(struct bgc2_header) == BGC2_HEADER_SIZE);
  input = check_fopen(filename, "rb");

  fread_fortran(hdr, BGC2_HEADER_SIZE, 1, input, 0);
  assert(hdr->magic == BGC_MAGIC);
  assert(hdr->version == 2);
  assert(hdr->format_group_data == GDATA_FORMAT_RMPVMAX);

  new_group_size = sizeof(GROUP_DATA_RMPVMAX)*((*num_groups)+hdr->ngroups);
  *groups = check_realloc(*groups, new_group_size, "Allocating groups.");
  fread_fortran((*groups) + (*num_groups), sizeof(GROUP_DATA_RMPVMAX), 
		hdr->ngroups, input, 0);
  *num_groups += hdr->ngroups;
  
  new_part_size = sizeof(PARTICLE_DATA_PV)*((*num_parts)+hdr->npart);
  *pdata = check_realloc(*pdata, new_part_size, "Allocating particles");
  p_start = 0;
  for (i=0; i<hdr->ngroups; i++) {
    fread_fortran((*pdata) + p_start, sizeof(PARTICLE_DATA_PV), 
		  groups[0][i].npart, input, 0);
    p_start += groups[0][i].npart;
  }
  *num_parts += hdr->npart;
  fclose(input);
}


