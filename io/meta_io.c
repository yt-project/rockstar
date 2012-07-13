#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../universal_constants.h"
#include "io_art.h"
#include "io_ascii.h"
#include "io_bgc2.h"
#include "io_gadget.h"
#include "io_generic.h"
#include "io_internal.h"
#include "io_tipsy.h"
#include "meta_io.h"
#include "../distance.h"
#include "../version.h"

char **snapnames = NULL;
char **blocknames = NULL;

void read_input_names(char *filename, char ***stringnames, int64_t *num_names) {
  int64_t i=0, j;
  char buffer[1024];
  FILE *input;
  char **names = NULL;

  if (!strlen(filename)) return;
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    while (strlen(buffer) && buffer[strlen(buffer)-1]=='\n')
      buffer[strlen(buffer)-1] = 0;
    if (!strlen(buffer)) continue;
    if (!(i%10)) names = check_realloc(names, sizeof(char *)*(i+10), 
					   "Allocating snapshot names.");
    names[i] = strdup(buffer);
    i++;
  }
  fclose(input);

  if (*stringnames) {
    for (j=0; j<*num_names; j++) free((*stringnames)[j]);
    free(*stringnames);
  }
  *num_names = i;
  *stringnames = names;
}

void get_input_filename(char *buffer, int maxlen, int64_t snap, int64_t block) {
  int64_t i=0, out=0, l=strlen(FILENAME);
  assert(snap < NUM_SNAPS);
  snprintf(buffer, maxlen, "%s/", INBASE);
  out=strlen(buffer);
  for (; (i<l)&&(out < (maxlen-1)); i++) {
    if (FILENAME[i] != '<') { buffer[out]=FILENAME[i]; buffer[out+1]=0; }
    else {
      if (!strncmp(FILENAME+i, "<snap>", 6)) {
	i+=5;
	if (snapnames) snprintf(buffer+out, maxlen-out, "%s", snapnames[snap]);
	else {
	  if (!strncasecmp(FILE_FORMAT, "GADGET", 6)) 
	    snprintf(buffer+out, maxlen-out, "%03"PRId64, snap);
	  else snprintf(buffer+out, maxlen-out, "%"PRId64, snap);
	}
      } 
      else if (!strncmp(FILENAME+i, "<block>", 7)) {
	i+=6;
	if (blocknames) snprintf(buffer+out, maxlen-out,"%s",blocknames[block]);
	else snprintf(buffer+out, maxlen-out, "%"PRId64, block);
      }
      else buffer[out] = FILENAME[i];
    }
    out = strlen(buffer);
  }
  buffer[out] = 0;
}

void get_output_filename(char *buffer, int maxlen, int64_t snap, int64_t chunk, char *type) {
  int64_t out = 0;
  snprintf(buffer, maxlen, "%s/", OUTBASE);
  out = strlen(buffer);
  if (snapnames) snprintf(buffer+out, maxlen-out, "halos_%s", snapnames[snap]);
  else snprintf(buffer+out, maxlen-out, "halos_%"PRId64, snap);
  out = strlen(buffer);
  snprintf(buffer+out, maxlen-out, ".%"PRId64".%s", chunk, type);
}


void read_particles(char *filename) {
  int64_t i, j, gadget = 0;
  int64_t p_start = num_p;
  float dx, ds, z, a, vel_mul;
  double *origin, origin_offset[3] = {0};
  if (!strcasecmp(FILE_FORMAT, "ASCII")) load_particles(filename, &p, &num_p);
  else if (!strncasecmp(FILE_FORMAT, "GADGET", 6)) {
    load_particles_gadget2(filename, &p, &num_p);
    gadget = 1;
  }
  else if (!strncasecmp(FILE_FORMAT, "ART", 3)) 
    load_particles_art(filename, &p, &num_p);
  else if (!strncasecmp(FILE_FORMAT, "INTERNAL", 8)) {
    load_particles_internal(filename, &p, &num_p);
  }
  else if (!strncasecmp(FILE_FORMAT, "GENERIC", 7)) {
    assert(load_particles_generic != NULL);
    load_particles_generic(filename, &p, &num_p);
  }
  else if (!strncasecmp(FILE_FORMAT, "TIPSY", 5)) {
    load_particles_tipsy(filename, &p, &num_p);
  }
  else {
    fprintf(stderr, "[Error] Unknown filetype %s!\n", FILE_FORMAT);
    exit(1);
  }

  if (LIMIT_RADIUS) {
    for (i=p_start; i<num_p; i++) {
      for (j=0, ds=0; j<3; j++) { dx = p[i].pos[j]-LIMIT_CENTER[j]; ds+=dx*dx; }
      if (ds > LIMIT_RADIUS*LIMIT_RADIUS) {
	num_p--;
	p[i] = p[num_p];
	i--;
      }
    }
  }

  if (LIGHTCONE) {
    init_cosmology();
    if (strlen(LIGHTCONE_ALT_SNAPS)) {
      for (i=0; i<3; i++)
	if (LIGHTCONE_ORIGIN[i] || LIGHTCONE_ALT_ORIGIN[i]) break;
      if (i<3) { //Same box coordinates, different intended locations
	if (LIGHTCONE == 1) {
	  for (i=0; i<3; i++) origin_offset[i] = LIGHTCONE_ORIGIN[i] - 
				LIGHTCONE_ALT_ORIGIN[i];
	}
      } else { //Offset everything
	for (i=0; i<3; i++) origin_offset[i] = -BOX_SIZE;
      }
      BOX_SIZE *= 2.0;
    }
    origin = (LIGHTCONE == 2) ? LIGHTCONE_ALT_ORIGIN : LIGHTCONE_ORIGIN;
    for (i=p_start; i<num_p; i++) {
      if (LIGHTCONE == 2) p[i].id = -p[i].id; //Make ids different
      for (j=0,dx=0; j<3; j++) {
	ds = p[i].pos[j] - origin[j];
	dx += ds*ds;
	p[i].pos[j] -= origin_offset[j];
      }
      if (!gadget) continue;
      dx = sqrt(dx);
      z = comoving_distance_h_to_redshift(dx);
      a = scale_factor(z);
      vel_mul = sqrt(a);
      for (j=0; j<3; j++) p[i].pos[j+3] *= vel_mul;
    }
  }
  output_config(NULL);
}

int _within_bounds(struct halo *h, float *bounds) {
  int64_t i;
  if (!bounds) return 1;
  for (i=0; i<3; i++) if (h->pos[i]<bounds[i]||h->pos[i]>bounds[i+3]) return 0;
  return 1;
}

int _should_print(struct halo *h, float *bounds) {
  if (!_within_bounds(h,bounds)) return 0;
  if ((h->num_p < MIN_HALO_OUTPUT_SIZE) ||
      (h->m * UNBOUND_THRESHOLD >= h->mgrav) ||
      ((h->m <= PARTICLE_MASS) && UNBOUND_THRESHOLD > 0)) return 0;
  return 1;
}

void print_ascii_header_info(FILE *output, float *bounds) {
  fprintf(output, "#a = %f\n", SCALE_NOW);
  if (bounds)
    fprintf(output, "#Bounds: (%f, %f, %f) - (%f, %f, %f)\n",
	    bounds[0], bounds[1], bounds[2], 
	    bounds[3], bounds[4], bounds[5]);
  fprintf(output, "#Om = %f; Ol = %f; h = %f\n", Om, Ol, h0);
  fprintf(output, "#Unbound Threshold: %f; FOF Refinement Threshold: %f\n",
	  UNBOUND_THRESHOLD, FOF_FRACTION);
  fprintf(output, "#Particle mass: %.5e Msun/h\n", PARTICLE_MASS);
  fprintf(output, "#Box size: %f Mpc/h\n", BOX_SIZE);
  fprintf(output, "#Total particles processed: %"PRId64"\n", num_p);
  fprintf(output, "#Force resolution assumed: %g Mpc/h\n", FORCE_RES);
  fprintf(output, "#Units: Masses in Msun / h\n"
	  "#Units: Positions in Mpc / h (comoving)\n"
	  "#Units: Velocities in km / s (physical, peculiar)\n"
	  "#Units: Radii in kpc / h (comoving)\n"
	  "#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n"
	  "#Units: Total energy in (Msun/h)*(km/s)^2 (physical)\n");
  fprintf(output, "#Rockstar Version: %s\n", ROCKSTAR_VERSION);
}

void output_ascii(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds) {
  char buffer[1024];
  int64_t i, id=0;
  FILE *output;
  struct halo *th;
  get_output_filename(buffer, 1024, snap, chunk, "ascii");
  //if (PARALLEL_IO) 
  output = check_fopen(buffer, "w");
  //else output = stdout;

  fprintf(output, "#id num_p m%s mbound_%s r%s vmax rvmax vrms x y z vx vy vz Jx Jy Jz E Spin PosUncertainty VelUncertainty bulk_vx bulk_vy bulk_vz BulkVelUnc n_core\n", MASS_DEFINITION, MASS_DEFINITION, MASS_DEFINITION);
  print_ascii_header_info(output, bounds);

  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    th = halos+i;
    fprintf(output, "%"PRId64" %"PRId64" %.3e %.3e"
	    " %f %f %f %f %f %f %f %f %f %f %g %g %g %g %g %f %f %f %f %f %f %"PRId64"\n", id+id_offset,
	    th->num_p, th->m, th->mgrav, th->r,	th->vmax, th->rvmax, th->vrms,
	    th->pos[0], th->pos[1], th->pos[2], th->pos[3], th->pos[4],
	    th->pos[5], th->J[0], th->J[1], th->J[2], th->energy, th->spin,
	    sqrt(th->min_pos_err), sqrt(th->min_vel_err), th->bulkvel[0],
	    th->bulkvel[1], th->bulkvel[2], sqrt(th->min_bulkvel_err),
	    th->n_core);
    id++;
  }
  //if (PARALLEL_IO) 
  fclose(output);
}

void print_child_particles(FILE *output, int64_t i, int64_t pid, int64_t eid) {
  int64_t j, child;
  struct particle *p2;
  for (j=0; j<halos[i].num_p; j++) {
    p2 = p + halos[i].p_start + j;
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %"PRId64" %"PRId64"\n", p2->pos[0], p2->pos[1], p2->pos[2], p2->pos[3], p2->pos[4], p2->pos[5], p2->id, i, pid, eid);
  }
  child = extra_info[i].child;
  while (child > -1) {
    if (halos[child].num_p < DOUBLE_COUNT_SUBHALO_MASS_RATIO*halos[pid].num_p)
      print_child_particles(output, child, pid, eid);
    child = extra_info[child].next_cochild;
  }
}

int compare_masses(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  if (c->m > d->m) return -1;
  if (d->m > c->m) return 1;
  return 0;
}

struct particle *output_p = NULL;
int64_t num_alloced_output_p = 0;

void add_more_output_p(int64_t total) {
  num_alloced_output_p = total+1000;
  output_p = check_realloc(output_p, sizeof(struct particle)*num_alloced_output_p,
			   "Allocating bound output particles.");
}

int64_t populate_copies(struct halo *h, float *cen, int64_t p_start)
{
  int64_t j, k, child;
  float dx = 0, r2, v2;
  struct particle *p2;

  if (num_alloced_output_p < p_start + h->num_p)
    add_more_output_p(p_start+h->num_p);
  memcpy(output_p + p_start, p + h->p_start, sizeof(struct particle)*h->num_p);
  for (j=0; j<h->num_p; j++) {
    output_p[p_start+j].id = h->p_start + j;
  }

  for (j=0; j<h->num_p; j++) {
    p2 = p + h->p_start + j;
    for (r2=0,k=0; k<3; k++) { dx = p2->pos[k]-cen[k]; r2+=dx*dx; }
    for (v2=0,k=3; k<6; k++) { dx = p2->pos[k]-cen[k]; v2+=dx*dx; }
    output_p[p_start + j].pos[0] = r2;
    output_p[p_start + j].pos[1] = v2;
  }
  child = extra_info[h-halos].child;
  p_start += h->num_p;
  while (child > -1) {
    p_start += populate_copies(halos+child, cen, p_start);
    child = extra_info[child].next_cochild;
  }
  return p_start;
}

int sort_radii(const void *a, const void *b) {
  const struct particle *c = a;
  const struct particle *d = b;
  if (c->pos[0] < d->pos[0]) return -1;
  if (d->pos[0] < c->pos[0]) return 1;
  return 0;
}

void print_bound_particles(struct halo *h, FILE *output, int64_t total_p) {
  int64_t i, num_bound = 0;
  double phi = 0;
  double potential_constant = 2.0*PARTICLE_MASS*Gc;
  struct particle *p2;
  for (i=0; i<total_p; i++)
    output_p[i].pos[0] = SCALE_NOW*sqrt(output_p[i].pos[0]);
  for (i=total_p-1; i>=0; i--) {
    if (potential_constant*(i/output_p[i].pos[0] + phi) > output_p[i].pos[1]) {
      p2 = p + output_p[i].id;
      if (!FULL_PARTICLE_CHUNKS)
	fprintf(output, "%"PRId64"\n", p2->id);
      else
	fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64"\n", p2->pos[0], p2->pos[1], p2->pos[2], p2->pos[3], p2->pos[4], p2->pos[5], p2->id, h->id);
      phi += 1.0/output_p[i].pos[0];
      num_bound++;
    }
  }
}

void output_particles(char *filename)
{
  char buffer[1024];
  int64_t i, j, id=0;
  FILE *output;
  struct particle *p2;
  sprintf(buffer, "%s/%s.list", OUTBASE, filename);
  output = check_fopen(buffer, "w");
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, NULL)) continue;
    for (j=0; j<halos[i].num_p; j++) {
      p2 = p + halos[i].p_start + j;
      fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64"\n", p2->pos[0], p2->pos[1], p2->pos[2], p2->pos[3], p2->pos[4], p2->pos[5], p2->id, id);
    }
    id++;
  }
  fclose(output);
}

void output_hmad(char *filename) {
  char buffer[1024];
  int64_t i, id=0;
  FILE *output;
  int64_t num_hp;

  qsort(halos, num_halos, sizeof(struct halo), compare_masses);
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, NULL)) continue;
    num_hp = populate_copies(halos+i,halos[i].pos,0);
    qsort(output_p, num_hp, sizeof(struct particle), sort_radii);

    sprintf(buffer, "%s/%s_%"PRId64".list", OUTBASE, filename, id);
    output = check_fopen(buffer, "w");

    if (!FULL_PARTICLE_CHUNKS)
      fprintf(output, "%f %f %f\n", halos[i].pos[0], halos[i].pos[1],
	      halos[i].pos[2]);
    print_bound_particles(halos+i, output, num_hp);
    id++;
    fclose(output);
  }
}

void output_full_particles(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds) {
  char buffer[1024];
  FILE *output;
  int64_t i, id=0;
  struct halo *th;

  if (chunk>=FULL_PARTICLE_CHUNKS) return;
  get_output_filename(buffer, 1024, snap, chunk, "particles");
  output = check_fopen(buffer, "w");

  fprintf(output, "#Halo table:\n");
  fprintf(output, "#id internal_id num_p m%s mbound_%s r%s vmax rvmax vrms x y z vx vy vz Jx Jy Jz energy spin\n", MASS_DEFINITION, MASS_DEFINITION, MASS_DEFINITION);
  fprintf(output, "#Particle table:\n");
  fprintf(output, "#x y z vx vy vz particle_id assigned_internal_haloid internal_haloid external_haloid\n");

  fprintf(output, "#Notes: As not all halos are printed, some halos may not have external halo ids.  (Hence the need to print internal halo ids).  Each particle is assigned to a unique halo; however, some properties (such as halo bound mass) are calculated including all substructure.  As such, particles belonging to subhalos are included in outputs; to exclude substructure, verify that the internal halo id is the same as the assigned internal halo id.\n");
  fprintf(output, "#Note also that halos with centers in particle ghost zones (outside the nominal boundaries) will not be printed, as their particle lists may be incomplete.\n");

  print_ascii_header_info(output, bounds);
  fprintf(output, "#Halo table begins here:\n");
  
  for (i=0; i<num_halos; i++) {
    th = halos+i;
    if (_should_print(th, bounds)) {
      th->id = id+id_offset;
      id++;
    } else { th->id = -1; }

    fprintf(output, "#%"PRId64" %"PRId64" %"PRId64" %.3e %.3e"
	    " %f %f %f %f %f %f %f %f %f %f %g %g %g %g %g\n", th->id, i,
	    th->num_p, th->m, th->mgrav, th->r, th->vmax, th->rvmax, th->vrms,
	    th->pos[0], th->pos[1], th->pos[2], th->pos[3], th->pos[4],
	    th->pos[5], th->J[0], th->J[1], th->J[2], th->energy, th->spin);
  }

  fprintf(output, "#Particle table begins here:\n");
  for (i=0; i<num_halos; i++) print_child_particles(output, i, i, halos[i].id);
  fclose(output);
  get_output_filename(buffer, 1024, snap, chunk, "particles");
  //gzip_file(buffer);
}

void delete_binary(int64_t snap, int64_t chunk) {
  char buffer[1024];
  get_output_filename(buffer, 1024, snap, chunk, "bin");
  unlink(buffer);
}

int64_t count_halos_to_print(float *bounds) {
  int64_t to_print = 0, i;
  for (i=0; i<num_halos; i++)
    if (_should_print(halos+i, bounds)) to_print++;
  return to_print;
}

void output_and_free_halos(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds) {
  if (!strcasecmp(OUTPUT_FORMAT, "BOTH") || !strcasecmp(OUTPUT_FORMAT, "ASCII"))
    output_ascii(id_offset, snap, chunk, bounds);
  if (!strcasecmp(OUTPUT_FORMAT, "BOTH")|| !strcasecmp(OUTPUT_FORMAT, "BINARY"))
    output_binary(id_offset, snap, chunk, bounds);

  if (chunk<FULL_PARTICLE_CHUNKS) {
    if (!fork()) {
      output_full_particles(id_offset, snap, chunk, bounds);
      exit(0);
    }
  }

  if (DUMP_PARTICLES[0] && (chunk >= DUMP_PARTICLES[1] &&
			    chunk <= DUMP_PARTICLES[2]))
    output_particles_internal(snap, chunk);

  output_bgc2(id_offset, snap, chunk, bounds);

  halos = check_realloc(halos, 0, "Freeing halo memory.");
  num_halos = 0;
}


void output_merger_catalog(int64_t snap, int64_t chunk, struct halo *halos, int64_t num_halos) {
  char buffer[1024];
  FILE *output;
  int64_t i, j;
  double m;
  struct halo *th;
  struct flock fl = {0};

  snprintf(buffer, 1024, "%s/out_%"PRId64".list", OUTBASE, snap);
  fl.l_type = F_WRLCK;
  fl.l_whence = SEEK_SET;
  fl.l_start = 0;
  fl.l_len = 0;

  if (chunk == 0) {
    output = check_fopen(buffer, "w");
    fcntl(fileno(output), F_SETLKW, &fl);
    fprintf(output, "#ID DescID M%s Vmax Vrms R%s Rs Np X Y Z VX VY VZ JX JY JZ Spin\n",
	    MASS_DEFINITION, MASS_DEFINITION);
    fprintf(output, "#a = %f\n", SCALE_NOW);
    fprintf(output, "#Om = %f; Ol = %f; h = %f\n", Om, Ol, h0);
    fprintf(output, "#Unbound Threshold: %f; FOF Refinement Threshold: %f\n",
	    UNBOUND_THRESHOLD, FOF_FRACTION);
    fprintf(output, "#Particle mass: %.5e Msun/h\n", PARTICLE_MASS);
    fprintf(output, "#Box size: %f Mpc/h\n", BOX_SIZE);
    fprintf(output, "#Units: Masses in Msun / h\n"
	    "#Units: Positions in Mpc / h (comoving)\n"
	    "#Units: Velocities in km / s (physical)\n"
	    "#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n"
	    "#Units: Radii in kpc / h (comoving)\n");
    fprintf(output, "#Rockstar Version: %s\n", ROCKSTAR_VERSION);
  }
  else {
    output = check_fopen(buffer, "a");
    fcntl(fileno(output), F_SETLKW, &fl);
  }
  
  for (i=0; i<num_halos; i++) {
    th = halos+i;
    if (LIGHTCONE) for (j=0; j<3; j++) th->pos[j] -= LIGHTCONE_ORIGIN[j];
    m = (BOUND_PROPS) ? th->mgrav : th->m;
    fprintf(output, "%"PRId64" %"PRId64" %.4e %.2f %.2f %.3f %.3f %"PRId64" %.5f %.5f %.5f %.2f %.2f %.2f %.3e %.3e %.3e %.5f\n",
	    th->id, th->desc, m, th->vmax, th->vrms, th->r, th->rs,
	    th->num_p, th->pos[0], th->pos[1], th->pos[2], th->pos[3],
	    th->pos[4], th->pos[5], th->J[0], th->J[1], th->J[2], th->spin);
  }
  fflush(output);
  fclose(output); //Also unlocks the file
}
