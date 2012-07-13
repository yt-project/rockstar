#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "bgc2.h"
#include "meta_io.h"
#include "io_util.h"
#include "../config_vars.h"
#include "../check_syscalls.h"
#include "../universal_constants.h"
#include "../rockstar.h"
#include "../groupies.h"

char **bgc2_snapnames = NULL;
int64_t num_bgc2_snaps = 0;
GROUP_DATA_RMPVMAX *gd = NULL;
float *bgc_group_center;

extern float particle_thresh_dens;

struct extended_particle {
  int64_t id, hid;
  float pos[6];
};

void populate_header(struct bgc2_header *hdr, int64_t id_offset, 
		     int64_t snap, int64_t chunk, float *bounds) {
  memset(hdr, 0, sizeof(struct bgc2_header));
  hdr->magic = BGC_MAGIC;
  hdr->version = 2;
  hdr->num_files = (PARALLEL_IO) ? NUM_WRITERS : 1;
  hdr->file_id = chunk;
  hdr->snapshot = snap;
  hdr->group_type = GTYPE_SO;
  hdr->format_part_data = PDATA_FORMAT_PV;
  hdr->format_group_data = GDATA_FORMAT_RMPVMAX;
  hdr->ngroups = 0;
  hdr->ngroups_total = 0;
  hdr->min_group_part = MIN_HALO_OUTPUT_SIZE;

  hdr->npart = 0;
  hdr->npart_total = 0;
  hdr->npart_orig = num_p;
  hdr->valid_part_ids = (!IGNORE_PARTICLE_IDS) ? 1 : 0;

  hdr->max_npart = 0;
  hdr->max_npart_total = 0;

  hdr->linkinglength = FOF_LINKING_LENGTH;
  calc_mass_definition();
  hdr->overdensity = particle_thresh_dens * PARTICLE_MASS / (Om * CRITICAL_DENSITY);
  hdr->time = SCALE_NOW;
  hdr->redshift = (SCALE_NOW>0) ? (1.0/(SCALE_NOW) - 1.0) : 1e10;
  hdr->Omega0 = Om;
  hdr->OmegaLambda = Ol;
  hdr->box_size = BOX_SIZE;

  if (bounds) for (int64_t i=0; i<6; i++) hdr->bounds[i] = bounds[i];

  hdr->Hubble0 = h0;
  hdr->GravConst = Gc;
  hdr->part_mass = PARTICLE_MASS;
}

void convert_to_extended_particles(struct extended_particle *ep) {
  int64_t i,j;
  p = (void *)ep;
  for (i=num_p-1; i>=0; i--) {
    ep[i].id = p[i].id;
    memmove(ep[i].pos, p[i].pos, sizeof(float)*6);
    ep[i].hid = -1;
  }
  for (i=0; i<num_halos; i++)
    for (j=0; j<halos[i].num_p; j++) ep[halos[i].p_start + j].hid = i;
}

float square_dist_from_center(struct extended_particle *a) {
  float ds=0, dx=0;
  int64_t i;
  for (i=0; i<3; i++) {
    dx = a->pos[i]-bgc_group_center[i]; ds += dx*dx;
  }
  return ds;
}

int sort_results_by_distance(const void *a, const void *b) {
  float c = square_dist_from_center(*((struct extended_particle **)a));
  float d = square_dist_from_center(*((struct extended_particle **)b));
  if (c < d) return -1;
  if (c > d) return 1;
  return 0;
}


#define FAST3TREE_TYPE struct extended_particle
#define FAST3TREE_PREFIX BGC2
#include "../fast3tree.c"

void output_bgc2(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds)
{
  char *buffer = NULL;
  int64_t i, j, k, id=0;
  FILE *output;
  struct bgc2_header *hdr = NULL;
  struct extended_particle *ep = NULL;
  struct fast3tree *ep_tree;
  struct fast3tree_results *ep_res;
  PARTICLE_DATA_PV *pd = NULL;
  int64_t num_to_print = count_halos_to_print(bounds);
  double dens_thresh;

  
  if (bgc2_snapnames == NULL) {
    if (!strlen(BGC2_SNAPNAMES)) return;
    read_input_names(BGC2_SNAPNAMES, &bgc2_snapnames, &num_bgc2_snaps);
  }

  for (i=0; i<num_bgc2_snaps; i++)
    if ((snapnames && !strcmp(snapnames[snap], bgc2_snapnames[i])) ||
	(!snapnames && atoi(bgc2_snapnames[i])==snap)) break;
  if (i==num_bgc2_snaps) return;

  
  assert(BGC2_HEADER_SIZE == sizeof(struct bgc2_header));
  buffer = check_realloc(buffer, 1025, "Allocating output buffer");
  get_output_filename(buffer, 1024, snap, chunk, "bgc2");
  output = check_fopen(buffer, "w");

  hdr = check_realloc(hdr, sizeof(struct bgc2_header),"Allocating BGC2 header");
  populate_header(hdr, id_offset, snap, chunk, bounds);
  hdr->ngroups = num_to_print;
  fwrite_fortran(hdr, BGC2_HEADER_SIZE, 1, output);
  if (!num_to_print) {
    fclose(output);
    return;
  }

  gd = check_realloc(gd, sizeof(GROUP_DATA_RMPVMAX)*num_to_print,
		     "Allocating output halo buffer");

  fwrite_fortran(gd, sizeof(GROUP_DATA_RMPVMAX), num_to_print, output);

  ep = check_realloc(p, sizeof(struct extended_particle)*num_p, "Allocating extended particle memory.");
  convert_to_extended_particles(ep);
  ep_tree = fast3tree_init(num_p, ep);
  ep_res = fast3tree_results_init();
  dens_thresh = particle_thresh_dens*(4.0*M_PI/3.0);

  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    
    fast3tree_find_sphere(ep_tree, ep_res, halos[i].pos, halos[i].r*1.1e-3);
    bgc_group_center = halos[i].pos;
    qsort(ep_res->points, ep_res->num_points, sizeof(struct extended_particle *), sort_results_by_distance);

    for (j=ep_res->num_points-1; j>=0; j--) {
      float r = sqrt(square_dist_from_center(ep_res->points[j]));
      if (r < FORCE_RES) r = FORCE_RES;
      float cur_dens = ((double)(j+1)/(r*r*r));
      if (cur_dens > dens_thresh) break;
    }
    if (j<0) j=0; //Almost certainly due to FP imprecision

    gd[id].id = id+id_offset;
    gd[id].parent_id = -1;
    gd[id].npart = j+1;
    gd[id].radius = cbrt((3.0/(4.0*M_PI))*(j+1)/particle_thresh_dens);
    gd[id].mass = (j+1)*PARTICLE_MASS;
    gd[id].vmax = halos[i].vmax;
    gd[id].rvmax = halos[i].rvmax;
    for (j=0; j<3; j++) {
      gd[id].pos[j] = halos[i].pos[j];
      gd[id].vel[j] = halos[i].pos[j+3];
    }
    gd[id].npart_self = 0;

    for (j=gd[id].npart-1; j>=(int64_t)gd[id].npart_self; j--) {
      struct extended_particle *tmp;
      if (ep_res->points[j]->hid == i) {
	tmp = ep_res->points[gd[id].npart_self];
	ep_res->points[gd[id].npart_self] = ep_res->points[j];
	ep_res->points[j] = tmp;
	gd[id].npart_self++;
	j++;
      }
    }

    hdr->npart += gd[id].npart;
    if (gd[id].npart > hdr->max_npart) {
      pd = check_realloc(pd, sizeof(struct particle)*gd[id].npart,
			 "Allocating particle output buffer.");
      hdr->max_npart = gd[id].npart;
    }

    for (j=0; j<gd[id].npart; j++) {
      pd[j].part_id = ep_res->points[j]->id;
      for (k=0; k<3; k++) {
	pd[j].pos[k] = ep_res->points[j]->pos[k];
	pd[j].vel[k] = ep_res->points[j]->pos[k+3];
      }
    }

    fwrite_fortran(pd, sizeof(PARTICLE_DATA_PV), gd[id].npart, output);
    id++;
  }
  rewind(output);
  fwrite_fortran(hdr, BGC2_HEADER_SIZE, 1, output);
  fwrite_fortran(gd, sizeof(GROUP_DATA_RMPVMAX), num_to_print, output);
  fclose(output);

  free(pd);
  free(buffer);
  free(hdr);
  gd = check_realloc(gd, 0, "Freeing group data.");
  fast3tree_results_free(ep_res);
  fast3tree_free(&ep_tree);
}


void load_bgc2_groups(char *filename, struct bgc2_header *hdr,
		    GROUP_DATA_RMPVMAX **groups, int64_t *num_groups)
{
  FILE *input;
  int64_t new_group_size;

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
  fclose(input);
}

