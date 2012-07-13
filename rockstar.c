/* The Rockstar Halo Finder.
   Copyright (C) 2011,2012  Peter Behroozi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   If so, it should be in a file called "LICENSE".
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include "rockstar.h"
#include "particle.h"
#include "fof.h"
#include "groupies.h"
#include "check_syscalls.h"
#include "config_vars.h"
#include "io/meta_io.h"

#define FAST3TREE_TYPE struct particle
#define FAST3TREE_PREFIX ROCKSTAR
#define POINTS_PER_LEAF 20
#include "fast3tree.c"
#define SKIP_THRESH 16

struct particle *p = NULL; 
int64_t num_p = 0;
struct fast3tree *tree = NULL;
struct fast3tree_results *rockstar_res = NULL;
char *skip = NULL;
struct fof *all_fofs = NULL;
int64_t num_all_fofs = 0;
int64_t num_fofs_tosend = 0, num_fofs_recvd = 0;

void rockstar(float *bounds, int64_t manual_subs) {
  int64_t i,j;
  float r;

  calc_mass_definition();
  r = AVG_PARTICLE_SPACING * FOF_LINKING_LENGTH;
  build_particle_tree();
  init_particle_smallfofs(num_p, p);
  skip = check_realloc(skip, sizeof(char)*num_p, "Allocating skip memory.");
  memset(skip, 0, sizeof(char)*num_p);
  
  for (i=0; i<num_p; i++) {
    if (skip[i]) continue;
    fast3tree_find_sphere(tree, rockstar_res, p[i].pos, r);
    link_particle_to_fof(p+i, rockstar_res->num_points, rockstar_res->points);
    if (rockstar_res->num_points > SKIP_THRESH) {
      for (j=0; j<rockstar_res->num_points; j++)
	skip[(rockstar_res->points[j]-p)] = 1;
      fast3tree_find_sphere(tree, rockstar_res, p[i].pos, 2.0*r);
      link_fof_to_fof(p+i, rockstar_res->num_points, rockstar_res->points);
    }
  }
  clear_particle_tree();
  skip = check_realloc(skip, 0, "Freeing skip memory.");
  build_fullfofs();
  all_fofs = return_fullfofs(&num_all_fofs);

  if (bounds) prune_fofs(bounds);

  if (!manual_subs) {
    for (i=0; i<num_all_fofs; i++)
      find_subs(all_fofs + i);

    rockstar_cleanup();
  } else {
    qsort(all_fofs, num_all_fofs, sizeof(struct fof), sort_fofs);
    num_fofs_tosend = num_all_fofs;
    num_fofs_recvd = 0;
  }
}

void find_unfinished_workunit(struct workunit_info *w, struct fof **fofs, struct particle **parts)
{
  int64_t i, np=0;
  struct fof *tf = *fofs;
  w->num_particles = w->num_fofs = w->num_halos = 0;
  w->seed = rand();
  for (; num_fofs_tosend>0 && (w->num_particles < LARGE_FOF); num_fofs_tosend--) {
    if (!(w->num_fofs % 1000))
      tf = *fofs = check_realloc(*fofs, sizeof(struct fof)*(w->num_fofs+1000),
				 "Allocating workunit fofs.");
    tf[w->num_fofs] = all_fofs[num_fofs_tosend-1];
    w->num_particles += tf[w->num_fofs].num_p;
    w->num_fofs++;
  }

  if (w->num_fofs == 1) {
    *parts = check_realloc(*parts, 0, "Freeing workunit particles.\n");
    *parts = tf[0].particles;
  } else {
    *parts = check_realloc(*parts, sizeof(struct particle)*w->num_particles,
			 "Allocating workunit particles.");
    for (i=0; i<w->num_fofs; i++) {
      memcpy(*parts + np, tf[i].particles, sizeof(struct particle)*tf[i].num_p);
      np+=tf[i].num_p;
    }
  }
}

void integrate_finished_workunit(struct workunit_info *w, struct fof *fofs, struct halo *h,
				 struct extra_halo_info *ei, struct particle *parts) {
  int64_t i, j=0, offset, np=0, hstart = num_halos;
  for (i=0; i<w->num_fofs; i++) {
    assert((fofs[i].particles >= p) && (fofs[i].particles < (p+num_p)));
    memcpy(fofs[i].particles, parts + np, sizeof(struct particle)*fofs[i].num_p);
    offset = (fofs[i].particles - p) - np;
    np += fofs[i].num_p;
    for (; (j < w->num_halos) && (h[j].p_start < np); j++) {
      add_new_halo();
      h[j].p_start += offset;
      if (ei[j].child >= 0) ei[j].child += hstart;
      if (ei[j].next_cochild >= 0) ei[j].next_cochild += hstart;
      if (ei[j].prev_cochild >= 0) ei[j].prev_cochild += hstart;
      if (ei[j].sub_of >= 0) ei[j].sub_of += hstart;
      halos[num_halos-1] = h[j];
      extra_info[num_halos-1] = ei[j];
    }
  }
  num_fofs_recvd += w->num_fofs;
}

int work_finished(void) {
  if (num_fofs_recvd >= num_all_fofs) return 1;
  return 0;
}

int sort_fofs(const void *a, const void *b) {
  const struct fof *c = a;
  const struct fof *d = b;
  if (c->num_p > LARGE_FOF) {
    if (c->num_p < d->num_p) return -1;
    if (c->num_p > d->num_p) return 1;
  }
  else if (d->num_p > LARGE_FOF) return -1;
  else {
    if (c->particles > d->particles) return -1;
    if (d->particles < c->particles) return 1;
  }
  return 0;
}

void do_workunit(struct workunit_info *w, struct fof *fofs) {
  int64_t i, processed_parts = 0;
  struct fof tmp;
  halos = check_realloc(halos, 0, "Freeing halo memory.");
  num_halos = 0;
  srand(w->seed);
  for (i=0; i<w->num_fofs; i++) {
    tmp = fofs[i];
    tmp.particles = p + processed_parts;
    find_subs(&tmp);
    processed_parts += tmp.num_p;
  }
}

void build_particle_tree(void) {
  int64_t i, dup_particles = 0;
  struct particle *last_p;
  tree = fast3tree_init(num_p, p);
  rockstar_res = fast3tree_results_init();
  if (IGNORE_PARTICLE_IDS)
    for (i=0; i<num_p; i++) p[i].id = i;

  if (num_p<2) return;
  last_p = p+(num_p-1);
  for (i=num_p-2; i>=0; i--) {
    if (!memcmp(p[i].pos, last_p->pos, sizeof(float)*6) || 
	last_p->id == p[i].id) {
      num_p--;
      p[i] = p[num_p];
      dup_particles++;
      continue;
    }
    last_p = p+i;
  }
  if (dup_particles) {
    fast3tree_rebuild(tree, num_p, p);
    if (dup_particles > 0.0001*num_p)
      fprintf(stderr, "[Warning] %"PRId64" duplicate particles removed.\n",
	      dup_particles);
  }
}

void clear_particle_tree(void) {
  fast3tree_results_free(rockstar_res);
  rockstar_res = NULL;
  fast3tree_free(&tree);
}

void rockstar_cleanup() {
  if (all_fofs) free(all_fofs);
  free_particle_copies();
  all_fofs = NULL;
  num_all_fofs = 0;
}

void prune_fofs(float *bounds) {
  int64_t i, j, k;
  for (i=0; i<num_all_fofs; i++) {
    for (j=0; j<all_fofs[i].num_p; j++) {
      for (k=0; k<3; k++)
	if (all_fofs[i].particles[j].pos[k] >= bounds[k] &&
	    all_fofs[i].particles[j].pos[k] <= bounds[k+3])
	  break;
      if (k<3) break;
    }
    if (j==all_fofs[i].num_p) {
      num_all_fofs--;
      all_fofs[i] = all_fofs[num_all_fofs];
      i--;
    }
  }
}

struct particle ** find_halo_sphere(struct halo *h, int64_t *num_results) {
  fast3tree_find_sphere(tree, rockstar_res, h->pos, h->r*1.1e-3);
  *num_results = rockstar_res->num_points;
  return rockstar_res->points;
}
