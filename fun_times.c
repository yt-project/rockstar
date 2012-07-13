#include <stdio.h>
#include <stdlib.h>
#include "io/meta_io.h"
#include "halo.h"
#include "particle.h"
#include "config_vars.h"
#include "check_syscalls.h"
#include "universe_time.h"
#include "litehash.h"
#include "rockstar.h"
#include "bounds.h"
#include "fun_times.h"

#define FAST3TREE_TYPE struct previous_halo
#define FAST3TREE_PREFIX FUN_TIMES
#include "fast3tree.c"

struct previous_halo *ph = NULL;
int64_t prev_snap = 0;
int64_t num_prev_halos = 0;
struct halo *prev_halo_buffer = NULL;
struct fast3tree *phtree = NULL;
struct fast3tree_results *phtree_results = NULL;
int64_t prev_pids[MAX_CORE_PARTICLES];
int64_t *hid_cache = NULL;

inline void add_to_previous_halos(struct halo *h, struct binary_output_header *bh) {
  if (!(num_prev_halos%1000))
    ph = check_realloc(ph, sizeof(struct previous_halo)*(num_prev_halos+1000),
		       "Allocating room for previous halos.");
  memcpy(ph[num_prev_halos].pos, h->pos, sizeof(float)*6);
  ph[num_prev_halos].m = h->m;
  ph[num_prev_halos].r = h->r;
  ph[num_prev_halos].num_p = h->num_p;
  ph[num_prev_halos].chunk = bh->chunk;
  ph[num_prev_halos].p_offset = h->p_start;
  ph[num_prev_halos].file_offset = sizeof(struct binary_output_header) +
    sizeof(struct halo)*bh->num_halos + h->p_start*sizeof(int64_t);
  num_prev_halos++;
}

void load_prev_binary_halos(int64_t snap, int64_t chunk, float *bounds, int64_t our_chunk) {
  FILE *input;
  char buffer[1024];
  struct binary_output_header bh;
  struct halo h;
  int64_t remaining = 0, to_read, i,j, h_start = num_prev_halos;
  double v_to_dx;

  get_output_filename(buffer, 1024, snap, chunk, "bin");
  input = check_fopen(buffer, "rb");
  check_fread(&bh, sizeof(struct binary_output_header), 1, input);
  assert(bh.magic == ROCKSTAR_MAGIC);
  assert(bh.num_halos >= 0);
  assert(bh.num_particles >= 0);

  //Conversion in Comoving Mpc/h / (km/s)
  //Note that the time units are in 1/H = 1/(h*100 km/s/Mpc)
  v_to_dx = 0.01*(scale_to_time(SCALE_NOW) - scale_to_time(bh.scale)) / 
    (0.5*(SCALE_NOW + bh.scale));

  remaining = bh.num_halos;
  while (remaining > 0) {
    to_read = PREV_HALO_BUFFER_SIZE;
    if (to_read > remaining) to_read = remaining;
    check_fread(prev_halo_buffer, sizeof(struct halo), to_read, input);
    remaining -= to_read;
    for (i=0; i<to_read; i++) {
      for (j=0; j<3; j++) 
	prev_halo_buffer[i].pos[j] += v_to_dx*prev_halo_buffer[i].pos[j+3];
      h = prev_halo_buffer[i];
      if (!bounds || _check_bounds(prev_halo_buffer[i].pos, h.pos, bounds) ||
	  our_chunk == chunk)
	add_to_previous_halos(&h, &bh);
    }
  }
  if (chunk == our_chunk) {
    hid_cache = check_realloc(hid_cache, sizeof(int64_t)*bh.num_particles,
			      "Allocating halo id cache");
    check_fread(hid_cache, sizeof(int64_t), bh.num_particles, input);
    for (i=h_start; i<num_prev_halos; i++)
      ph[i].file_offset = -1;
  }
  fclose(input);
}

void load_previous_halos(int64_t snap, int64_t chunk, float *bounds) {
  int64_t rchunk;
  struct binary_output_header bh;
  float overlap_region[6];

  num_prev_halos = 0;
  prev_snap = snap-1;
  if (!phtree) {
    phtree = fast3tree_init(num_prev_halos, ph);
    phtree_results = fast3tree_results_init();
  }
  else fast3tree_rebuild(phtree, num_prev_halos, ph);
  if (snap == STARTING_SNAP || LIGHTCONE || !PARALLEL_IO) return;
  prev_halo_buffer = check_realloc(prev_halo_buffer,
				   sizeof(struct halo)*PREV_HALO_BUFFER_SIZE,
				   "Allocating previous halo buffer.");
  for (rchunk = 0; rchunk < NUM_WRITERS; rchunk++) {
    load_binary_header(snap-1, rchunk, &bh);
    if (!bounds || bounds_overlap(bh.bounds, bounds, overlap_region, OVERLAP_LENGTH))
      load_prev_binary_halos(snap-1, rchunk, bounds, chunk);
  }
  prev_halo_buffer = check_realloc(prev_halo_buffer, 0,
				   "Freeing previous halo buffer.");
  fast3tree_rebuild(phtree, num_prev_halos, ph);
}

int64_t prev_halo_acceptable(struct halo *h, struct previous_halo *tph) {
  int64_t j;
  float ds, dx=0, dv=0;
  if (tph->num_p < h->num_p*0.5) return 0;
  for (j=0; j<6; j++) {
    ds = h->pos[j]-tph->pos[j];
    if (j < 3) dx += ds*ds;
    else dv += ds*ds;
  }
  if (dv > h->vrms*h->vrms) return 0;
  return 1;
}

int sort_previous_halos(const void *a, const void *b) {
  struct previous_halo *c = *((struct previous_halo **)a);
  struct previous_halo *d = *((struct previous_halo **)b);
  if (c->chunk < d->chunk) return -1;
  if (c->chunk > d->chunk) return 1;
  if (c->file_offset < d->file_offset) return -1;
  if (c->file_offset > d->file_offset) return 1;
  return 0;
}

int sort_core_particles(const void *a, const void *b) {
  const struct particle *c = a;
  const struct particle *d = b;
  if (c->pos[0] < d->pos[0]) return -1;
  if (c->pos[0] > d->pos[0]) return 1;
  return 0;
}

void convert_and_sort_core_particles(struct halo *h, struct particle *hp) {
  int64_t i, j;
  float ds, dx;
  for (i=0; i<h->num_p; i++) {
    for (j=0,dx=0; j<3; j++) { ds = h->pos[j]-hp[i].pos[j]; dx+= ds*ds; }
    hp[i].pos[0] = ds;
  }
  qsort(hp, h->num_p, sizeof(struct particle), sort_core_particles);
  for (i=0; i<h->num_p; i++)
    hp[i].pos[0] = p[hp[i].id].pos[0];
}


float find_previous_mass(struct halo *h, struct particle *hp, int64_t *best_num_p) {
  int64_t i, j, last_chunk = -1, max_particles, best_particles = 0;
  int64_t cur_part, remaining;
  char buffer[1024];
  struct previous_halo *tph, *best_ph=NULL;
  struct litehash *lh = NULL;
  FILE *input = NULL;

  *best_num_p = 0;
  if (h->num_p < 100 || !num_prev_halos) return 0;
  fast3tree_find_sphere(phtree, phtree_results, h->pos, h->r);
  if (!phtree_results->num_points) return 0;
  for (i=0; i<phtree_results->num_points; i++) {
    if (!prev_halo_acceptable(h, phtree_results->points[i])) {
      phtree_results->num_points--;
      phtree_results->points[i] = phtree_results->points[phtree_results->num_points];
      i--;
    }
  }
  if (!phtree_results->num_points) return 0;

  qsort(phtree_results->points, phtree_results->num_points, sizeof(struct previous_halo *),
	sort_previous_halos);

  max_particles = MAX_CORE_PARTICLES;
  if (max_particles >= h->num_p) max_particles = h->num_p;
  else convert_and_sort_core_particles(h, hp);
  
  lh = new_litehash(8);
  for (i=0; i<max_particles; i++)
    lh_setval(lh, &(p[hp[i].id].id), (void *)1);

  for (i=0; i<phtree_results->num_points; i++) {
    tph = phtree_results->points[i];
    remaining = tph->num_p;
    cur_part = 0;
    if (remaining > MAX_CORE_PARTICLES) remaining = MAX_CORE_PARTICLES;

    if (tph->file_offset >= 0) {
      if (tph->chunk != last_chunk) {
	if (last_chunk > -1) fclose(input);
	get_output_filename(buffer, 1024, prev_snap, tph->chunk, "bin");
	input = check_fopen(buffer, "rb");
	last_chunk = tph->chunk;
      }
      check_fseeko(input, tph->file_offset, SEEK_SET);
      check_fread(prev_pids, sizeof(int64_t), remaining, input);
    } else {
      memcpy(prev_pids, hid_cache+tph->p_offset, sizeof(int64_t)*remaining);
    }

    for (j=0; j<remaining; j++)
      if (lh_getval(lh, prev_pids + j)) cur_part++;

    if (cur_part > best_particles) {
      best_particles = cur_part;
      best_ph = tph;
    }
  }

  if (last_chunk > -1) fclose(input);
  free_litehash(lh);
  if (best_ph && best_ph->m > 1e13 && best_particles > max_particles*0.1) {
    fprintf(stderr, "Hnow: %f %f %f (%"PRId64"; %e); Phalo: %f %f %f (%e); %"PRId64"\n",
	    h->pos[0], h->pos[1], h->pos[2], h->num_p, h->num_p*PARTICLE_MASS,
	    best_ph->pos[0], best_ph->pos[1], best_ph->pos[2], best_ph->m,
	    best_particles);
  }
  *best_num_p = best_particles;
  if (best_ph && (best_particles > max_particles*0.1)) return best_ph->m;
  return 0;
}
