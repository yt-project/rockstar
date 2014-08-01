#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "../config_vars.h"
#include "../config.h"
#include "../check_syscalls.h"
#include "../rockstar.h"
#include "../version.h"
#include "../io/meta_io.h"
#include "../io/bgc2.h"
#include "../io/io_bgc2.h"
#include "SDF.h"

extern GROUP_DATA_RMPVMAX *gd;

#define RMPVMAX_DESC \
"/* Positions in Mpc/h, Velocities km/s, masses Msol/h */\n\
struct {\n\
    int64_t id;\n\
    int64_t parent_id;\n\
    int64_t npart;\n\
    int64_t npart_self; /* Excluding substructure */\n\
    float radius;\n\
    float mass;\n\
    float x;\n\
    float y;\n\
    float z;\n\
    float vx;\n\
    float vy;\n\
    float vz;\n\
    float vmax;\n\
    float rvmax;\n\
}"

static int rev_sort_by_mass(const void *a, const void *b) {
  const GROUP_DATA_RMPVMAX *c = a;
  const GROUP_DATA_RMPVMAX *d = b;
  if (c->mass > d->mass) return -1;
  if (c->mass < d->mass) return 1;
  return 0;
}

int main(int argc, char **argv)
{
    int64_t i, snap=-1, did_config = 0;
    struct bgc2_header *hdrs = NULL;
    int64_t gnobj;
    float pmass;
    double vol0_hinvMpc;
    char buffer[1024];
    char filename[256];
    FILE *output;
    int bin, nbins;
    int ilogm_start = 8;
    int ilogm_end = 17;
    int bins_per_decade = 10;
    int64_t *counts;
    double *bin_center, *bin_left, *bin_right;

    if (argc < 2) {
	printf("Usage: %s [-c config] [-s snapnum] [-b blocknum]\n", argv[0]);
	exit(1);
    }
    
    for (i=1; i<argc-1; i++) {
	if (!strcmp("-c", argv[i])) { do_config(argv[i+1]); i++; did_config=1; }
	if (!strcmp("-s", argv[i])) { snap = atoi(argv[i+1]); i++; }
    }
    if (!did_config) do_config(NULL);
    if (strlen(SNAPSHOT_NAMES)) 
	read_input_names(SNAPSHOT_NAMES, &snapnames, &NUM_SNAPS);
    if (strlen(BLOCK_NAMES))
	read_input_names(BLOCK_NAMES, &blocknames, &NUM_BLOCKS);
    
    hdrs = (struct bgc2_header *)calc_bgc2_parents(snap);
    
    gnobj = hdrs[0].ngroups_total;
    vol0_hinvMpc = hdrs[0].box_size*hdrs[0].box_size*hdrs[0].box_size;
    pmass = hdrs[0].part_mass;

    /* Remove those below minimum (keep subhalos) */
    for (i = 0; i < gnobj; i++) {
	if (gd[i].mass/pmass < MIN_HALO_OUTPUT_SIZE-0.001) {
	    gd[i--] = gd[--gnobj];
	}
    }
    if (gnobj == 0) exit(0);
    gd = check_realloc(gd, gnobj*sizeof(GROUP_DATA_RMPVMAX), "Realloc group data.");
    /* Sort from largest to smallest */
    qsort(gd, gnobj, sizeof(GROUP_DATA_RMPVMAX), rev_sort_by_mass);

    strncpy(filename, FILENAME, sizeof(filename));
    for (i = 0; i < strlen(filename); i++) {
	if (filename[i] == '<') filename[i] = '\0';
    }
    snprintf(buffer, 1024, "%s/%s%s.top1000", OUTBASE, filename, snapnames[snap]);
    output = check_fopen(buffer, "w");
    fprintf(output, "#ID DescID M%s Vmax Vrms R%s Rs Np X Y Z VX VY VZ Parent_ID\n",
	    MASS_DEFINITION2, MASS_DEFINITION2);
    fprintf(output, "# overdensity=%.3f redshift=%.3f part_mass=%g\n",
	    hdrs[0].overdensity, hdrs[0].redshift, hdrs[0].part_mass);
    for (i = 0; i < ((gnobj > 1000) ? 1000 : gnobj); i++) {
	fprintf(output, "%"PRId64" -1 %.6g %.3f 0 %.3f 0 %"PRId64" %f %f %f %f %f %f %"PRId64"\n", gd[i].id, gd[i].mass, gd[i].vmax, gd[i].radius*1.0e3, gd[i].npart, gd[i].pos[0],  gd[i].pos[1],  gd[i].pos[2], gd[i].vel[0], gd[i].vel[1], gd[i].vel[2], gd[i].parent_id);
    }
    fclose(output);

    snprintf(buffer, 1024, "%s/%sbgc2_%s", OUTBASE, filename, snapnames[snap]);
    SDFwrite(buffer, gnobj, gnobj, gd, sizeof(GROUP_DATA_RMPVMAX), RMPVMAX_DESC,
	     "nhalos", SDF_INT64, gnobj, 
	     "R0", SDF_DOUBLE, BOX_SIZE/2.0,
	     "L0", SDF_DOUBLE, BOX_SIZE,
	     "BOX_SIZE", SDF_DOUBLE, BOX_SIZE,
	     "a", SDF_DOUBLE, SCALE_NOW,
	     "SCALE_NOW", SDF_DOUBLE, SCALE_NOW,
	     "overdensity", SDF_DOUBLE, hdrs[0].overdensity,
	     "redshift", SDF_DOUBLE, 1.0/SCALE_NOW-1.0,
	     "part_mass", SDF_DOUBLE, hdrs[0].part_mass,
	     "Omega0_m", SDF_DOUBLE, Om,
	     "Omega0_lambda", SDF_DOUBLE, Ol,
	     "h_100", SDF_DOUBLE, h0,
	     "m200b", SDF_INT, 1, 
	     "units_rockstar", SDF_INT, 1, 
	     "Rockstar_version", SDF_STRING, Rockstar_version,
	     "libSDF_version", SDF_STRING, libSDF_version,
	       NULL);

    nbins = (ilogm_end-ilogm_start)*bins_per_decade;
    counts = check_realloc(0, nbins*sizeof(int64_t), "alloc bin counts.");
    bin_center = check_realloc(0, nbins*sizeof(double), "alloc bin center.");
    bin_left = check_realloc(0, nbins*sizeof(double), "alloc bin left.");
    bin_right = check_realloc(0, nbins*sizeof(double), "alloc bin right.");
    for (i = 0; i < nbins; i++) {
	double x;
	bin = i+ilogm_start*bins_per_decade;
	/* place bins between particle masses */
	x = floor(pow(10.0, bin/(double)bins_per_decade-0.5/bins_per_decade)/pmass);
	bin_left[i] = pow(10.0, (log10(x*pmass)+log10((x+1.0)*pmass))/2.0);
	x = floor(pow(10.0, bin/(double)bins_per_decade+0.5/bins_per_decade)/pmass);
	bin_right[i] = pow(10.0, (log10(x*pmass)+log10((x+1.0)*pmass))/2.0);
	bin_center[i] = pow(10.0, (log10(bin_left[i])+log10(bin_right[i]))/2.0);
	counts[i] = 0;
    }
    
    /* loop assumes gd sorted by mass */
    for (i = 0, bin = nbins-1; i < gnobj; i++) {
	double m = gd[i].mass;
	if (gd[i].parent_id != -1) continue;
	if (bin_right[bin] < m) {
	    fprintf(stderr, "mass outside of maximum bin %d %g %g\n", bin, m, bin_right[bin]);
	    exit(1);
	}
	while (m < bin_left[bin]) {
	    bin--;
	    if (!bin) {
		fprintf(stderr, "mass outside of minimum bin %g %g\n", m, bin_left[bin]);
		exit(1);
	    }
	}
	if (m > bin_right[bin]) {
	    fprintf(stderr, "bins inconsistent %g %g\n", bin_left[bin], bin_right[bin]);
	    exit(1);
	}
	counts[bin]++;
    }

    snprintf(buffer, 1024, "%s/%s%s.hist_m200b", OUTBASE, filename, snapnames[snap]);
    output = check_fopen(buffer, "w");

    fprintf(output, "# bin_center dn/dlog10M n bin_left bin_right\n");
    fprintf(output, "# volume %.8g [Mpc^3/h^3]\n", vol0_hinvMpc);
    fprintf(output, "# bins_per_decade %d\n", bins_per_decade);
    fprintf(output, "# overdensity %f\n", hdrs[0].overdensity);
    fprintf(output, "# minimum_npart %ld\n", MIN_HALO_OUTPUT_SIZE);
    for (i = 0; i < nbins; i++) {
	if (counts[i] && bin_left[i]/pmass >= MIN_HALO_OUTPUT_SIZE) {
	    fprintf(output, "%14.8lg %14.8lg %9ld %14.8lg %14.8lg\n", 
		    bin_center[i], counts[i]*bins_per_decade/vol0_hinvMpc, counts[i], bin_left[i], bin_right[i]);
	}
    }
    fclose(output);
    exit(0);
}
