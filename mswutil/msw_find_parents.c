#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../check_syscalls.h"
#include "../io/stringparse.h"
#include "../universal_constants.h"
#include "../version.h"
#include <math.h>

#define EXTRA_HALO_INFO int64_t descid, np; \
  float alt_m[4], J[3], spin, bullock_spin, Xoff, Voff, \
    b_to_a, c_to_a, A[3], klypin_rs, kin_to_pot, m_all;
#include "read_tree.h"

double BOX_SIZE=250;
double SCALE_NOW, Om, Ol, h0, PARTICLE_MASS;

struct halo_list all_halos = {0};
struct halo_tree halo_tree = {0};

#define GROUP_LIST all_halos.halos
#define RADIUS r
#define FAST3TREE_TYPE struct halo
#include "../fast3tree.c"
#define parent pid
#include "../parents.c"
#undef parent

double vir_density(double a) {
  double x = 1.0/(1.0+a*a*a*Ol/Om)-1.0;
  return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x));
}

float calc_dens_definition(char *md) {
  int64_t length = strlen(md);
  char last_char = (length) ? md[length-1] : 0;
  float matter_fraction = 1.0/(1.0+pow(SCALE_NOW, 3)*Ol/Om);
  float cons = Om * CRITICAL_DENSITY;
  char *mass = md;
  float thresh_dens;
  if (mass[0] == 'm' || mass[0] == 'M') mass++;

  if (last_char == 'b' || last_char == 'B')
    thresh_dens = atof(mass) * cons;
  else if (last_char == 'c' || last_char == 'C')
    thresh_dens = atof(mass) * cons / matter_fraction;
  else {
    if (strcasecmp(md, "vir")) md = "vir";
    thresh_dens = vir_density(SCALE_NOW) * cons;
  }
  return thresh_dens;
}

void read_hlist(char *filename) {
  int64_t n;
  FILE *input;
  struct xhalo h = {0};
  char buffer[1024];

  SHORT_PARSETYPE;
  #define NUM_INPUTS 33
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F, F, F,  //  #id desc_id mvir vmax vrms
      F, F, D64, F,       //  Rvir Rs Np x 
      F, F, F, F, F,      // y z vx vy vz 
      F, F, F, F, F,      // JX JY JZ Spin rs_klypin
      F, F, F, F, F,      // M_all M1 M2 M3 M4
      F, F, F, F, F,      // Xoff Voff spin_bullock b_to_a c_to_a 
      F, F, F, F,         // A[x] A[y] A[z] T/|U|
    };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(h.id),
                            &(h.descid),
			    &(h.mvir), &(h.vmax), &(h.vrms), &(h.rvir), &(h.rs), 
			    &(h.np), 
			    &(h.pos[0]), &(h.pos[1]), &(h.pos[2]), 
			    &(h.vel[0]), &(h.vel[1]), &(h.vel[2]),
			    &(h.J[0]), &(h.J[1]), &(h.J[2]), &(h.spin),
			    &(h.klypin_rs), &(h.m_all), &(h.alt_m[0]), 
			    &(h.alt_m[1]), &(h.alt_m[2]), &(h.alt_m[3]),
			    &(h.Xoff), &(h.Voff), &(h.bullock_spin), 
			    &(h.b_to_a), &(h.c_to_a), &(h.A[0]), 
			    &(h.A[1]), &(h.A[2]), &(h.kin_to_pot)};

  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') {
	if (!strncmp("#Box size:", buffer, 9)) {
	    if (sscanf(buffer, "#Box size: %lf Mpc/h", &BOX_SIZE) != 1) {
		fprintf(stderr, "# failed to parse Box size\n");
	    } else {
		fprintf(stderr,  "# parsed Box size of %lf from input\n", BOX_SIZE);
	    }
	}
	if (!strncmp("#Particle mass:", buffer, 9)) {
	    if (sscanf(buffer, "#Particle mass: %lf Msun/h", &PARTICLE_MASS) != 1) {
		fprintf(stderr, "# failed to parse Particle mass\n");
	    } else {
		fprintf(stderr,  "# parsed Particle mass of %lf from input\n", PARTICLE_MASS);
	    }
	}
	if (!strncmp("#a =", buffer, 4)) {
	    if (sscanf(buffer, "#a = %lf", &SCALE_NOW) != 1) {
		fprintf(stderr, "# failed to parse a\n");
	    } else {
		fprintf(stderr,  "# parsed a of %lf from input\n", SCALE_NOW);
	    }
	}
	if (!strncmp("#Om =", buffer, 4)) {
	    if (sscanf(buffer, "#Om = %lf; Ol = %lf; h = %lf;", &Om, &Ol, &h0) != 3) {
		fprintf(stderr, "# failed to parse Om, Ol, h\n");
	    } else {
		fprintf(stderr,  "# parsed Om = %lf; Ol = %lf; h0 = %lf from input\n", Om, Ol, h0);
	    }
	}
    }
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    if (n<NUM_INPUTS) continue;
    if (!(all_halos.num_halos%300000))
      all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*(all_halos.num_halos+300000), "Allocating Halos.");

    for (int i=0; i<3; i++) {
	all_halos.halos[all_halos.num_halos].pos[i] = h.pos[i];
	all_halos.halos[all_halos.num_halos].vel[i] = h.vel[i];
    }
    all_halos.halos[all_halos.num_halos].r = h.rvir;
    all_halos.halos[all_halos.num_halos].m200b = h.alt_m[0];
    all_halos.halos[all_halos.num_halos].m2500c = h.alt_m[3];
    all_halos.halos[all_halos.num_halos].vmax = h.vmax;
    all_halos.halos[all_halos.num_halos].spin = h.spin;
    all_halos.halos[all_halos.num_halos].kin_to_pot = h.kin_to_pot;
    all_halos.halos[all_halos.num_halos].id = h.id;
    all_halos.halos[all_halos.num_halos].pid = h.pid;
    all_halos.num_halos++;
  }
  fclose(input);
  
  all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*all_halos.num_halos, "Allocating Halos.");

  char MASS_DEFINITION[] = "m200b";
  double thresh_dens = calc_dens_definition(MASS_DEFINITION)*(4.0*M_PI/3.0);
  /* switch to r200b */
  for (n=0; n<all_halos.num_halos; n++) {
    struct halo *th = all_halos.halos + n;
    th->r = cbrt(th->m200b/thresh_dens)*1e3;
  }

  find_parents(all_halos.num_halos);

  unsigned char bytes[4] = {0x12, 0x34, 0x56, 0x78};
  unsigned int *byteorder = (unsigned int *)&bytes[0];

  printf("# SDF-1.0\n");
  printf("parameter byteorder = 0x%x;\n", *byteorder);
  printf("int64_t nhalos = %ld;\n", all_halos.num_halos);
  printf("double BOX_SIZE = %lf;\n", BOX_SIZE);
  printf("double SCALE_NOW = %lf;\n", SCALE_NOW);
  printf("double overdensity = %lf;\n", 200.0);
  printf("double redshift = %lf;\n", 1.0/SCALE_NOW-1.0);
  printf("double part_mass = %lg;\n", PARTICLE_MASS);
  printf("double Omega0 = %lf;\n", Om);
  printf("double Lambda_prime = %lf;\n", Ol);
  printf("double h_100 = %lf;\n", h0);
  printf("int so200b = 1;\n");
  printf("int rockstar_units = 1;\n");
  printf("char Rockstar_version[] = \"%s\";\n", Rockstar_version);
  printf("struct {\n");
  printf("    float x, y, z;\n");
  printf("    float vx, vy, vz;\n");
  printf("    float m200b, m2500c, vmax, spin, kin_to_pot;\n");
  printf("    int64_t id, pid;\n");
  printf("}[%ld];\n", all_halos.num_halos);
  printf("\f\n# SDF-EOH\n");
  for (n=0; n<all_halos.num_halos; n++) {
      struct halo *th = all_halos.halos + n;
      fwrite(th->pos, 3, sizeof(float), stdout);
      fwrite(th->vel, 3, sizeof(float), stdout);
      fwrite(&th->m200b, 1, sizeof(float), stdout);
      fwrite(&th->m2500c, 1, sizeof(float), stdout);
      fwrite(&th->vmax, 1, sizeof(float), stdout);
      fwrite(&th->spin, 1, sizeof(float), stdout);
      fwrite(&th->kin_to_pot, 1, sizeof(float), stdout);
      fwrite(&th->id, 1, sizeof(int64_t), stdout);
      fwrite(&th->pid, 1, sizeof(int64_t), stdout);
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


