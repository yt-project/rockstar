#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include "rockstar.h"
#include "halo.h"
#include "fof.h"
#include "particle.h"
#include "groupies.h"
#include "subhalo_metric.h"
#include "check_syscalls.h"
#include "config_vars.h"
#include "universal_constants.h"
#include "potential.h"
#include "nfw.h"
#include "distance.h"
#include "fun_times.h"
#include "jacobi.h"

#define FAST3TREE_DIM 6
#define POINTS_PER_LEAF 50
#define FAST3TREE_PREFIX GROUPIES
#define FAST3TREE_TYPE struct particle
#include "fast3tree.c"

struct particle *copies = NULL; //For storing phase-space FOFs
int64_t *particle_halos = NULL;
float *particle_r = NULL;
struct potential *po = NULL;
int64_t num_alloc_pc = 0;

struct fof *subfofs = NULL;
int64_t num_subfofs = 0, num_alloced_subfofs = 0;

int64_t num_halos = 0;
struct halo *halos = NULL;
struct extra_halo_info *extra_info = NULL;

struct fast3tree_results *res = NULL;
struct fast3tree *phasetree = NULL;

int num_alloc_gh = 0;
struct halo **growing_halos = NULL;

int64_t *halo_ids = NULL;
int64_t num_alloced_halo_ids = 0;

float particle_thresh_dens = 0, particle_rvir_dens = 0;

double vir_density(double a) {
  double x = 1.0/(1.0+a*a*a*Ol/Om)-1.0;
  return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x));
}

void calc_mass_definition(void) {
  int64_t length = strlen(MASS_DEFINITION);
  char last_char = (length) ? MASS_DEFINITION[length-1] : 0;
  float matter_fraction = 1.0/(1.0+pow(SCALE_NOW, 3)*Ol/Om);
  float cons = Om * CRITICAL_DENSITY / PARTICLE_MASS; // background density
  char *mass = MASS_DEFINITION;
  if (mass[0] == 'm' || mass[0] == 'M') mass++;

  if (last_char == 'b' || last_char == 'B')
    particle_thresh_dens = atof(mass) * cons;
  else if (last_char == 'c' || last_char == 'C')
    particle_thresh_dens = atof(mass) * cons / matter_fraction;
  else {
    if (strcasecmp(MASS_DEFINITION, "vir")) MASS_DEFINITION = "vir";
    particle_thresh_dens = vir_density(SCALE_NOW) * cons;
  }
  particle_rvir_dens = vir_density(SCALE_NOW) * cons;
}

void lightcone_set_scale(float *pos) {
  int64_t i;
  float ds=0, dx=0, z;
  for (i=0; i<3; i++) { ds = pos[i]-LIGHTCONE_ORIGIN[i]; dx+=ds*ds; }
  z = comoving_distance_h_to_redshift(sqrt(dx));
  SCALE_NOW = scale_factor(z);
  calc_mass_definition();
}

void add_new_halo(void) {
  int i;
  if ((num_halos % 1000)==0) {
    halos = check_realloc(halos, sizeof(struct halo)*(num_halos+1000),
			  "Allocating room for halos.");
    extra_info = check_realloc(extra_info, sizeof(struct extra_halo_info)*(num_halos+1000),
			  "Allocating room for extra halo info.");
    memset(halos+num_halos, 0, sizeof(struct halo)*1000);
    for (i=num_halos; i<num_halos+1000; i++) {
      extra_info[i].child = extra_info[i].next_cochild = 
	extra_info[i].prev_cochild = extra_info[i].sub_of = -1;
      extra_info[i].max_metric = 0;
      halos[i].flags |= GROWING_FLAG;
    }
  }
  num_halos++;
}

void add_more_halo_ids(void) {
  num_alloced_halo_ids+=1000;
  halo_ids = check_realloc(halo_ids, sizeof(int64_t)*num_alloced_halo_ids,
			   "Allocating room for halo ids.");
}

void add_more_growing_halos(void) {
  num_alloc_gh += 1000;
  growing_halos = check_realloc(growing_halos, sizeof(struct halo *)
				*num_alloc_gh,
				"Allocating room for growing halos.");
}

void calc_basic_halo_props(struct halo *h) {
  int64_t j, k;
  double pos[6] = {0}, pos2[6] = {0}, x;
  double pos_err, vel_err;
  h->r = h->vrms = 0;
  for (j=0; j<h->num_p; j++)
    for (k=0; k<6; k++) pos[k] += copies[h->p_start + j].pos[k];

  for (k=0; k<6; k++) pos[k] /= (double)h->num_p;

  for (j=0; j<h->num_p; j++)
    for (k=0; k<6; k++) {
      x = copies[h->p_start + j].pos[k] - pos[k];
      pos2[k] += x*x;
    }

  for (k=0; k<6; k++) {
    if (k<3) h->r += pos2[k] / (double)h->num_p;
    else h->vrms += pos2[k] / (double)h->num_p;
  }

  pos_err = h->r / (double)h->num_p;
  vel_err = h->vrms / (double)h->num_p;

  if ((!h->min_pos_err) || (h->min_pos_err > pos_err)) {
    h->min_pos_err = pos_err;
    h->n_core = h->num_p;
    for (k=0; k<3; k++) h->pos[k] = pos[k];
  }

  if ((!h->min_vel_err) || (h->min_vel_err > vel_err)) {
    h->min_vel_err = vel_err;
    for (k=3; k<6; k++) h->pos[k] = pos[k];
  }
  for (k=3; k<6; k++) h->bulkvel[k-3] = pos[k];

  h->m = h->num_p;
  if (!h->num_child_particles) h->num_child_particles = h->num_p;
  
  h->r = cbrt(h->num_p/((4.0*M_PI/3.0)*particle_rvir_dens));
  h->child_r = cbrt(h->num_child_particles/((4.0*M_PI/3.0)*particle_rvir_dens));
  //h->r = h->child_r;
  h->vrms = sqrt(h->vrms);
}

void add_ang_mom(double L[3], float c[6], float pos[6]) {
  // L = r x p;
#define cross(a,x,y,s) L[a] s (pos[x]-c[x])*(pos[y+3]-c[y+3])
  cross(0,1,2,+=);
  cross(0,2,1,-=);
  cross(1,2,0,+=);
  cross(1,0,2,-=);
  cross(2,0,1,+=);
  cross(2,1,0,-=);
#undef cross
}

int dist_compare(const void *a, const void *b) {
  float c = ((struct potential *)a)->r2;
  float d = ((struct potential *)b)->r2;
  if (c>d) return 1;
  if (c<d) return -1;
  return 0;
}

void _reset_potentials(struct halo *base_h, struct halo *h, float *cen, int64_t p_start, int64_t level, int64_t potential_only) {
  int64_t j, k;
  float dx, r2;
  memset(po + p_start, 0, sizeof(struct potential)*h->num_p);
  for (j=0; j<h->num_p; j++) {
    r2 = 0;
    for (k=0; k<3; k++) { dx=copies[h->p_start+j].pos[k] - cen[k]; r2+=dx*dx; }
    po[p_start + j].r2 = r2;
    memcpy(po[p_start+j].pos, copies[h->p_start+j].pos, sizeof(float)*6);
    if (potential_only) po[p_start + j].ke = -1;
    if (h==base_h) po[p_start + j].flags = 1;
    if (!potential_only && (h->num_p < base_h->num_p*0.03))
      po[p_start+j].flags = 2;
  }
}

int64_t calc_particle_radii(struct halo *base_h, struct halo *h, float *cen, int64_t p_start, int64_t level, int64_t potential_only) {
  int64_t j, total_p = p_start, child, first_child, parent; 
  int64_t do_potential_only = 0;

  //Break accidental graph loops
  if (level >= num_alloced_halo_ids) add_more_halo_ids();
  halo_ids[level] = h-halos;
  for (j=0; j<level; j++) if (halo_ids[j] == halo_ids[level]) return p_start;

  _reset_potentials(base_h, h, cen, p_start, level, potential_only);

  first_child = child = extra_info[h-halos].child;
  total_p += h->num_p;
  while (child > -1) {
    if (!potential_only)
      do_potential_only = (halos[child].num_child_particles < 
          DOUBLE_COUNT_SUBHALO_MASS_RATIO*base_h->num_child_particles) ? 0 : 1;
    else do_potential_only = 1;
    total_p = calc_particle_radii(base_h, halos + child,
                                  cen, total_p, level+1, do_potential_only);
    child = extra_info[child].next_cochild;
    assert(child != first_child);
  }

  parent = extra_info[h-halos].sub_of;
  if ((h == base_h) && (parent > -1) &&
    (halos[parent].num_child_particles*INCLUDE_HOST_POTENTIAL_RATIO < h->num_child_particles)){
    total_p = calc_particle_radii(base_h, halos + parent,
                                  cen, total_p, level+1, 1);
  }
  return total_p;
}


/*int64_t calc_particle_radii2(struct halo *base_h, struct halo *h, int64_t p_start, int64_t potential_only) {
  struct halo_metric **children;
  struct halo *child, *par = NULL;
  int64_t num_children = 0, i, total_p = p_start, parent;
  _reset_potentials(h, h, h->pos, total_p, 0, potential_only);
  total_p += h->num_p;

  parent = extra_info[h-halos].sub_of;
  if (parent > -1) par = halos + parent;
  children = find_children(h, par, h->child_r*2, &num_children);
  for (i=0; i<num_children; i++) {
    child = children[i]->target;
    if ((child == base_h) || (potential_only &&
       (extra_info[child-halos].cur_parent == (base_h - halos)))) continue;
    _reset_potentials(h, child, h->pos, total_p, 0, potential_only);
    total_p += child->num_p;
    extra_info[child-halos].cur_parent = h - halos;
    }*/

  /*  if ((h == base_h) && (parent > -1) &&
    (halos[parent].num_child_particles*INCLUDE_HOST_POTENTIAL_RATIO < h->num_child_particles))
    total_p = calc_particle_radii2(base_h, halos + parent, total_p, 1);
  */
/*  return total_p;
    }*/


void _calc_num_child_particles(struct halo *h) {
  int64_t child, first_child;

  if (h->num_child_particles) return;
  h->num_child_particles = h->num_p;

  first_child = child = extra_info[h-halos].child;
  while (child > -1) {
    _calc_num_child_particles(halos + child);
    if (halos[child].num_p < DOUBLE_COUNT_SUBHALO_MASS_RATIO*h->num_p)
      h->num_child_particles += halos[child].num_child_particles;
    child = extra_info[child].next_cochild;
    assert(child != first_child);
  }
}

void calc_num_child_particles(int64_t h_start) {
  int64_t i;
  for (i=h_start; i<num_halos; i++) halos[i].num_child_particles = 0;
  for (i=h_start; i<num_halos; i++) 
    if (!halos[i].num_child_particles) _calc_num_child_particles(halos + i);
}

void calculate_corevel(struct halo *h, struct potential *po, int64_t total_p) {
  //Assumes po is already sorted.
  int64_t i, j;
  double vel[3]={0};
  int64_t core_max, velthresh, rvir_max;
  double var[3]={0}, thisvar, bestvar=0;
  double rvir_thresh = particle_rvir_dens*(4.0*M_PI/3.0);

  for (j=total_p-1; j>=0; j--)
    if (j*j > (po[j].r2*po[j].r2*po[j].r2)*(rvir_thresh*rvir_thresh)) break;
  rvir_max = j;
  for (j=total_p-1; j>=0; j--)
    if (j*j > (po[j].r2*po[j].r2*po[j].r2)*(100*rvir_thresh*rvir_thresh)) break;
  core_max = j;
  
  velthresh = pow(h->vrms / 5.0, 2.0);
  if (core_max < velthresh) core_max = velthresh;

  for (i=0; i<rvir_max; i++) {
    for (j=0; j<3; j++) {
      double delta = po[i].pos[j+3] - vel[j];
      vel[j] += delta / ((double)(i+1));
      var[j] += delta * (po[i].pos[j+3]-vel[j]);
    }
    thisvar = (var[0]+var[1]+var[2]);
    if ((i < 10) || (thisvar < bestvar*(i-3)*i)) {
      if (i > 3) bestvar = thisvar / (double)((i-3)*i);
      else bestvar = 0;
      if (i < core_max) {
	h->n_core = i;
	h->min_vel_err = bestvar;
	for (j=0; j<3; j++) h->corevel[j] = vel[j];
      }
      for (j=0; j<3; j++) h->bulkvel[j] = vel[j];
      h->min_bulkvel_err = bestvar;
    }
  }
}


void _calc_additional_halo_props(struct halo *h, int64_t total_p, int64_t bound)
{
  int64_t j, k, part_mdelta=0, num_part=0;
  double dens_thresh = particle_thresh_dens*(4.0*M_PI/3.0);
  double rvir_thresh = particle_rvir_dens*(4.0*M_PI/3.0);
  double vmax_conv = PARTICLE_MASS/SCALE_NOW;
  double r, circ_v, vmax=0, rvmax=0, energy = 0, L[3] = {0}, Jh, m=0, de;
  double vrms[3]={0}, vavg[3]={0}, mdiff;
  double cur_dens;

  for (j=0; j<total_p; j++) {
    if (bound && (po[j].pe < po[j].ke)) continue;
    num_part++;
    r = sqrt(po[j].r2);
    if (r < FORCE_RES) r = FORCE_RES;
    cur_dens = ((double)num_part/(r*r*r));
    
    if (cur_dens > dens_thresh) {
      part_mdelta = num_part;
      add_ang_mom(L, h->pos, po[j].pos);
      de = po[j].ke - po[j].pe;
      if (isfinite(de)) energy += de;
      for (k=0; k<3; k++) { //Calculate Vrms
	mdiff = po[j].pos[k+3]-vavg[k];
	vavg[k] += mdiff/(double)num_part;
	vrms[k] += mdiff*(po[j].pos[k+3]-vavg[k]);
      }
    }

    if (cur_dens > rvir_thresh) {
      circ_v = (double)num_part/r;
      if (part_mdelta && circ_v > vmax) {
	vmax = circ_v;
	rvmax = r;
      }
    }
  }

  m = part_mdelta*PARTICLE_MASS;
  if (!bound) h->m = m;
  else h->mgrav = m;
  for (k=0; k<3; k++) vrms[k] = (vrms[k] > 0) ? (vrms[k]/part_mdelta) : 0;
  if ((bound && BOUND_PROPS) || !(bound || BOUND_PROPS)) {
    h->vrms = sqrt(vrms[0] + vrms[1] + vrms[2]); 
    h->vmax = VMAX_CONST*sqrt(vmax*vmax_conv);
    h->rvmax = rvmax*1e3;
    h->r = cbrt((3.0/(4.0*M_PI))*part_mdelta/particle_thresh_dens)*1e3;
    h->rs = calc_scale_radius(h->m, h->r, h->vmax, h->rvmax, SCALE_NOW);
    for (j=0; j<3; j++) h->J[j] = PARTICLE_MASS*SCALE_NOW*L[j];
    h->energy = energy * PARTICLE_MASS * Gc / SCALE_NOW;
    Jh = PARTICLE_MASS*SCALE_NOW*sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
    h->spin = (m>0) ? (Jh * sqrt(fabs(h->energy)) / (Gc*pow(m, 2.5))) : 0;
  }
}


//Assumes center + velocity already calculated.
void calc_additional_halo_props(struct halo *h) {
  int64_t j, total_p;
  double dens_thresh;

  if (LIGHTCONE) lightcone_set_scale(h->pos);
  dens_thresh = particle_thresh_dens*(4.0*M_PI/3.0);
  if (h->num_p < 1) return;
  total_p = //calc_particle_radii2(h, h, 0, 0);
    calc_particle_radii(h, h, h->pos, 0, 0, 0);
  if (BOUND_OUT_TO_HALO_EDGE) {
    qsort(po, total_p, sizeof(struct potential), dist_compare);
    for (j=total_p-1; j>=0; j--)
      if (j*j / (po[j].r2*po[j].r2*po[j].r2) > dens_thresh*dens_thresh) break;
    if (total_p) total_p = j+1;
  }

  if (total_p>1) compute_potential(po, total_p);
  for (j=0; j<total_p; j++)
    if (po[j].ke < 0) {
      total_p--;
      po[j] = po[total_p];
      j--;
    }
  qsort(po, total_p, sizeof(struct potential), dist_compare);
  calculate_corevel(h, po, total_p);
  if (extra_info[h-halos].sub_of > -1)
    compute_kinetic_energy(po, total_p, h->corevel, h->pos);
  else
    compute_kinetic_energy(po, total_p, h->bulkvel, h->pos);

  _calc_additional_halo_props(h, total_p, 0);
  _calc_additional_halo_props(h, total_p, 1);
}

void _find_subfofs_at_r(struct fof *f, float target_r) {
  int64_t i;
  init_particle_smallfofs(f->num_p, f->particles);
  for (i=0; i<f->num_p; i++) {
    fast3tree_find_sphere(phasetree, res, f->particles[i].pos, target_r);
    link_particle_to_fof(f->particles + i, res->num_points, res->points);
  }
  build_fullfofs();
}

#define MAX_PARTICLES_TO_SAMPLE 10000
void _find_subfofs_better2(struct fof *f,  float thresh) {
  int64_t i, j, num_test = MAX_PARTICLES_TO_SAMPLE;
  float target_r = 0;
  norm_sd(f, thresh);
  fast3tree_rebuild(phasetree, f->num_p, f->particles);
  if (num_test > f->num_p) num_test = f->num_p;
  for (i=0; i<num_test; i++) {
    if (f->num_p <= MAX_PARTICLES_TO_SAMPLE) j = i;
    else { j = rand(); j<<=31; j+=rand(); j%=(f->num_p); }
    particle_r[i] = fast3tree_find_next_closest_distance(phasetree, res, 
			     f->particles[j].pos);
  }
  target_r = find_median_r(particle_r, num_test, thresh);
  _find_subfofs_at_r(f, target_r);
}

void reassign_halo_particles(int64_t p_start, int64_t p_end) {
  int64_t last_halo, j;
  partition_sort_particles(p_start, p_end, copies, particle_halos);
  last_halo = particle_halos[p_start];
  halos[last_halo].p_start = p_start;
  for (j=p_start + 1; j<p_end; j++) {
    if (particle_halos[j] == last_halo) continue;
    halos[last_halo].num_p = j - halos[last_halo].p_start;
    last_halo = particle_halos[j];
    halos[last_halo].p_start = j;
  }
  halos[last_halo].num_p = j - halos[last_halo].p_start;
}

int could_be_poisson_or_force_res(struct halo *h1, struct halo *h2, int64_t *is_force_res) {
  float dx, r=0, v=0, mpe, mve; //, vt1, vt2;
  int64_t k;
  *is_force_res = 0;
  if (!h1->min_pos_err || !h1->min_vel_err) return 1;
  for (k=0; k<3; k++) {
    dx = h1->pos[k]-h2->pos[k];      r+=dx*dx;
    dx = h1->pos[k+3]-h2->pos[k+3];  v+=dx*dx;
  }
  mpe = h1->min_pos_err;
  mve = h1->min_vel_err;
  dx = (r / mpe + v / mve) / 2.0;
  if (!(dx > 100)) return 1;
  
  r = sqrt(r);
  v = sqrt(v);
  if ((h1->r+h2->r > r) && (1.5*(h1->vrms+h2->vrms)) > v &&
      (r < 1.5*FORCE_RES)) {
    *is_force_res = 1;
    return 1;
  }
  return 0;
}

int64_t _find_biggest_parent(int64_t h_start, int use_temporal_info) {
  int64_t i, max_i = h_start, max_p, num_m1, num_m2;
  float m1 = -1;
  for (i=h_start+1; i<num_halos; i++)
    if (halos[i].r > halos[max_i].r)
      max_i = i;

  max_p = halos[max_i].num_p;
  if (use_temporal_info && TEMPORAL_HALO_FINDING && PARALLEL_IO) {
    for (i=h_start; i<num_halos; i++) {
      if (i==max_i) continue;
      if (max_p*0.25 < halos[i].num_p) {
	if (m1 < 0)
	  m1 = find_previous_mass(halos+max_i, copies+halos[max_i].p_start,
				  &num_m1);
	float m2 = find_previous_mass(halos+i, copies + halos[i].p_start,
				      &num_m2);
	if (m1 && m2 && ((m2 > m1) || ((m2 == m1) && (num_m2 > num_m1)))) {
	  max_i = i;
	  m1 = m2;
	  num_m1 = num_m2;
	}
      }
    }
  }  
  return max_i;
}

void _fix_parents(int64_t h_start) {
  int64_t i, sub_of, num_m1, num_m2;
  float m1, m2;

  for (i=h_start; i<num_halos; i++) {
    extra_info[i].next_cochild = extra_info[i].prev_cochild = 
      extra_info[i].child = -1;
    //extra_info[i].prev_mass = 0;
    if (extra_info[i].sub_of == i) extra_info[i].sub_of = -1;
  }

  if (TEMPORAL_HALO_FINDING && PARALLEL_IO) {
    for (i=h_start; i<num_halos; i++) {
      sub_of = extra_info[i].sub_of;
      if (sub_of == i) sub_of = extra_info[i].sub_of = -1;
      if (sub_of > -1 && halos[i].num_p > 0.25*halos[sub_of].num_p) {
	m2 = find_previous_mass(halos+i, copies+halos[i].p_start, &num_m2);
	m1 = find_previous_mass(halos+sub_of, copies+halos[sub_of].p_start,
				&num_m1);
	if (m1 && m2 && ((m2 > m1) || ((m2 == m1) && (num_m2 > num_m1)))) {
	  extra_info[i].sub_of = extra_info[sub_of].sub_of;
	  extra_info[sub_of].sub_of = i;
	  if (extra_info[sub_of].sub_of == sub_of) 
	    extra_info[sub_of].sub_of = -1;
	  i--;
	}
      }
    }    
  }

  for (i=h_start; i<num_halos; i++) {
    extra_info[i].max_metric = 1e10;
    sub_of = extra_info[i].sub_of;
    if (sub_of > -1)
      extra_info[i].max_metric = _calc_halo_dist(halos+i, halos + sub_of);      
    else continue;
    int64_t next_child = extra_info[sub_of].child;
    extra_info[i].next_cochild = next_child;
    extra_info[i].prev_cochild = -1;
    extra_info[sub_of].child = i;
    if (next_child > -1)
      extra_info[next_child].prev_cochild = i;
  }
}

void output_level(int64_t p_start, int64_t p_end, int64_t h_start, int64_t level)
{
  int64_t i;
  char buffer[1024];
  snprintf(buffer, 1024, "%s/levels", OUTBASE);
  FILE *output = check_fopen(buffer, "a");
  for (i=p_start; i<p_end; i++) {
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %"PRId64"\n",
	    copies[i].pos[0], copies[i].pos[1], copies[i].pos[2], 
	    copies[i].pos[3], copies[i].pos[4], copies[i].pos[5], 
	    p[copies[i].id].id, particle_halos[i], level);
  }
  fclose(output);

  snprintf(buffer, 1024, "%s/halos.levels", OUTBASE);
  output = check_fopen(buffer, "a");
  for (i=h_start; i<num_halos; i++) {
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %f %f %f %f %"PRId64" %"PRId64" %"PRId64"\n",
	    halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], 
	    halos[i].bulkvel[0], halos[i].bulkvel[1], halos[i].bulkvel[2], 
	    halos[i].num_p, halos[i].num_child_particles, 
	    halos[i].r, halos[i].vrms, sqrt(halos[i].min_pos_err), 
	    sqrt(halos[i].min_vel_err), i, extra_info[i].sub_of, level);
  }
  fclose(output);
}

void _find_subs(struct fof *f, int level) {
  int64_t f_start, f_end, h_start, i, j, f_index;
  int64_t p_start, num_growing_halos = 0, max_i = 0, is_force_res;

  //Find subFOFs
  p_start = f->particles - copies;
  f_index = f - subfofs;  
  _find_subfofs_better2(f, FOF_FRACTION);
  f_start = num_subfofs;
  copy_fullfofs(&subfofs, &num_subfofs, &num_alloced_subfofs);
  f_end = num_subfofs;

  h_start = num_halos;
  for (i=f_start; i<f_end; i++)
    if (subfofs[i].num_p > MIN_HALO_PARTICLES)
      _find_subs(subfofs + i, level+1);

  //Convert particle positions back to normal
  if (level>0) f = subfofs + f_index;
  for (j=0; j<f->num_p; j++) particle_halos[p_start + j] = -1;
  for (i=h_start; i<num_halos; i++)
    for (j=0; j<halos[i].num_p; j++)
      particle_halos[halos[i].p_start + j] = i;

  for (j=0; j<f->num_p; j++) {
    struct particle *c = f->particles + j;
    memcpy(c->pos, p[c->id].pos, sizeof(float)*6);
  }

  if (h_start == num_halos) add_new_halo(); //New seed halo
  max_i = _find_biggest_parent(h_start, 1);

  num_growing_halos=1;
  if (num_growing_halos >= num_alloc_gh) add_more_growing_halos();
  halos[max_i].flags |= GROWING_FLAG;
  extra_info[max_i].sub_of = -1;
  growing_halos[0] = halos + max_i;
  for (i=h_start; i<num_halos; i++) {
    if ((i==max_i) || !(halos[i].flags & GROWING_FLAG)) continue;
    if (could_be_poisson_or_force_res(halos+i, halos+max_i, &is_force_res)) {
      halos[i].flags -= (halos[i].flags & GROWING_FLAG);
      extra_info[i].sub_of = max_i;
      if (!is_force_res) {
	for (j=0; j<halos[i].num_p; j++)
	  particle_halos[halos[i].p_start+j] = max_i;
	halos[i].num_p = 0;
      }
      continue;
    }
    if (num_growing_halos >= num_alloc_gh) add_more_growing_halos();
    growing_halos[num_growing_halos] = halos+i;
    num_growing_halos++;
  }

  if (num_growing_halos==1) {
    for (j=0; j<f->num_p; j++)
      if (particle_halos[p_start + j] < 0)
	particle_halos[p_start + j] = max_i;
  } else {
    build_subtree(growing_halos, num_growing_halos);
    for (j=0; j<f->num_p; j++) {
      if (particle_halos[p_start + j] < 0) {
	struct halo *h = find_best_halo(copies+p_start+j, halos+max_i);
	particle_halos[p_start + j] = h - halos;
	while (extra_info[h-halos].sub_of > -1) {
	  float max_metric = extra_info[h-halos].max_metric;
	  if (max_metric < 3) max_metric = 3;
	  if (calc_particle_dist(h, copies+p_start+j) > max_metric) {
	    particle_halos[p_start + j] = extra_info[h-halos].sub_of;
	    h = halos + extra_info[h-halos].sub_of;
	  }
	  else break;
	}
      }
    }
  }

  reassign_halo_particles(p_start, p_start + f->num_p);
  //calc_num_child_particles(h_start);
  for (i=0; i<num_growing_halos; i++) calc_basic_halo_props(growing_halos[i]);
  max_i = _find_biggest_parent(h_start, 0);
  for (i=0; i<num_growing_halos; i++)
    extra_info[growing_halos[i]-halos].sub_of = 
      find_best_parent(growing_halos[i], halos+max_i) - halos;
  _fix_parents(h_start);
  calc_num_child_particles(h_start);
  for (i=0; i<num_growing_halos; i++) calc_basic_halo_props(growing_halos[i]);

  if (OUTPUT_LEVELS) output_level(p_start, p_start+f->num_p, h_start, level);


  num_subfofs = f_start;
  if (!level) {
    num_growing_halos = num_halos - h_start;
    while (num_alloc_gh < num_growing_halos) add_more_growing_halos();
    for (i=0; i<num_growing_halos; i++) growing_halos[i] = halos+h_start+i;
    for (i=0; i<num_growing_halos; i++) {
      calc_basic_halo_props(growing_halos[i]);
      convert_and_sort_core_particles(growing_halos[i], 
				      copies + growing_halos[i]->p_start);
    }
    build_subtree(growing_halos, num_growing_halos);
    max_i = _find_biggest_parent(h_start, 0);
    for (i=0; i<num_growing_halos; i++)
      extra_info[growing_halos[i]-halos].sub_of = 
	find_best_parent(growing_halos[i], halos+max_i) - halos;
    _fix_parents(h_start);
  }
}

void find_subs(struct fof *f) {
  struct fof cf;
  int64_t i, h_start = num_halos;

  if (!phasetree) phasetree = fast3tree_init(0, NULL);
  if (!res) res = fast3tree_results_init();

  if (f->num_p > num_alloc_pc) alloc_particle_copies(f->num_p);
  if (!f->num_p) return;
  memcpy(copies, f->particles, sizeof(struct particle)*f->num_p);
  for (i=0; i<f->num_p; i++) copies[i].id = (f->particles-p)+i; //Hijack particle IDs
  cf = *f;
  cf.particles = copies;

  if (LIGHTCONE) lightcone_set_scale(f->particles->pos);
  num_subfofs = 0;
  _find_subs(&cf, 0);
  num_subfofs = 0;
  for (i=0; i<f->num_p; i++) copies[i] = p[copies[i].id];
  calc_num_child_particles(h_start);
  for (i=h_start; i<num_halos; i++) calc_basic_halo_props(halos + i);
  for (i=h_start; i<num_halos; i++) calc_additional_halo_props(halos + i);

  memcpy(f->particles, copies, sizeof(struct particle)*f->num_p);
  for (i=h_start; i<num_halos; i++)
    halos[i].p_start += (f->particles - p);
}


void alloc_particle_copies(int64_t total_copies) {
  int64_t max_particle_r = MAX_PARTICLES_TO_SAMPLE;
  if (total_copies - num_alloc_pc < 1000) total_copies = num_alloc_pc + 1000;
  copies = check_realloc(copies, sizeof(struct particle)*total_copies,
			 "Allocating room for particle copies.");
  particle_halos = check_realloc(particle_halos, sizeof(int64_t)*total_copies,
			 "Allocating room for particle halo-links.");
  if (max_particle_r > total_copies) max_particle_r = total_copies;
  particle_r = check_realloc(particle_r, sizeof(float)*max_particle_r,
			 "Allocating room for particle radii.");
  po = check_realloc(po, sizeof(struct potential)*total_copies,
			"Allocating room for particle distances.");
  num_alloc_pc = total_copies;
}

void free_particle_copies(void) {
  num_alloc_pc = 0;
  copies = check_realloc(copies, 0, "Freeing copies.");
  particle_halos = check_realloc(particle_halos, 0, "Freeing particle links.");
  particle_r = check_realloc(particle_r, 0, "Freeing particle radii.");
  po = check_realloc(po, 0, "Freeing potentials.");
  free_subtree();
}

int64_t rad_partition(float *rad, int64_t left, int64_t right, int64_t pivot_ind) {
  float pivot = rad[pivot_ind], tmp;
  int64_t si, i;
#define SWAP(a,b) { tmp = rad[a]; rad[a] = rad[b]; rad[b] = tmp; }
  SWAP(pivot_ind, right-1);
  si = right-2;
  for (i = left; i<si; i++) {
    if (rad[i] > pivot) { SWAP(i, si); si--; i--; }
  }
  if (rad[si] < pivot) si++;
  SWAP(right-1, si);
  return si;
#undef SWAP
}

inline float random_unit(void) {
  return(((float)(rand()%(RAND_MAX))/(float)(RAND_MAX)));
}

float find_median_r(float *rad, int64_t num_p, float frac) {
  int64_t pivot_index, k = num_p * frac;
  int64_t left = 0, right = num_p;
  assert(num_p>0);
  if (num_p < 2) return rad[0];
  while (1) {
    pivot_index = rad_partition(rad, left, right, 
				left + random_unit()*(right-left));
    if (k == pivot_index || (rad[left]==rad[right-1])) return rad[k];
    if (k < pivot_index) right = pivot_index;
    else left = pivot_index+1;
  }
}


void norm_sd(struct fof *f, float thresh) {
  double pos[6] = {0};
  double corr[6][6];
  double sig_x, sig_v;
  int64_t i,j,k;

  if (!f->num_p) return;

  for (i=0; i<f->num_p; i++)
    for (j=0; j<6; j++) pos[j]+=f->particles[i].pos[j];

  for (j=0; j<6; j++) pos[j]/=(double)f->num_p;
  for (j=0; j<6; j++) for (k=0; k<6; k++) corr[j][k] = 0;

  for (i=0; i<f->num_p; i++)
    for (j=0; j<6; j++) {
      f->particles[i].pos[j] -= pos[j];
      for (k=j; k<6; k++)
	corr[j][k]+=f->particles[i].pos[j]*f->particles[i].pos[k];
    }

  for (j=0; j<6; j++)
    for (k=j; k<6; k++)
      corr[j][k]/=(double)f->num_p;

  calc_deviations(corr, &sig_x, &sig_v);

  if (!sig_x || !sig_v) return;

  for (i=0; i<f->num_p; i++)
    for (j=0; j<6; j++)
      f->particles[i].pos[j] /= ((j < 3) ? sig_x : sig_v);
}
