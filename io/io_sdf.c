/* Copyright (C) 2013-2014  Michael S. Warren
   Explicit permission granted to distribute with Rockstar under the GPLv3
   The SDF library is available at https://bitbucket.org/JohnSalmon/SDF
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "io_sdf.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../distance.h"
#include "../hubble.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"
#include "SDF.h"

#define one_kpc 3.08567802e16 /* km */
#define one_Gyr 3.1558149984e16 /* sec */

#define Error(...) { fprintf(stderr, __VA_ARGS__); exit(1); }
#define singlError(...) { if (block == 0) fprintf(stderr, __VA_ARGS__); exit(1); }

void sdf_split(int64_t gnobj, int nproc, int procnum, int *nobj, int64_t *start)
{
    int leftover;

    *nobj = gnobj / nproc;
    *start = (int64_t)(*nobj) * procnum;
    leftover = gnobj - (*nobj) * nproc;
    if (procnum < leftover) {
	++(*nobj);
	*start += procnum;
    } else {
	*start += leftover;
    }
}


void sdf_check(struct particle *p, int64_t p_start, int64_t nelems) {
    int64_t i, j;
    int allzero, nwarn = 0;

    for (i=0; i<nelems; i++) {
	allzero = 1;
	for (j=0; j<6; j++) {
	    if (p[p_start+i].pos[j] != 0.0) allzero = 0;
	}
	if (allzero){
	    fprintf(stderr, "[Warning] sdf_check particle %ld all zeros\n", 
		    p[p_start+i].id);
	    nwarn++;
	    if (nwarn > 10) Error("[Error] Bad input file\n");
	}
    }
}


void sdf_renumber_particles(struct particle *p, int64_t p_start, int64_t nelems, 
			    int64_t start) {
    int64_t i;
    fprintf(stderr, "[Warning] sdf_renumber %ld to %ld\n", 
	    p_start+start, p_start+start+nelems);
    for (i=0; i<nelems; i++) {
	p[p_start+i].id = start+p_start+i;
    }
}

void sdf_rescale_particles(struct particle *p, int64_t p_start, int64_t nelems, 
			   double length_conversion, double vel_conversion) {
  int64_t i, j;
  double pos_rescale = length_conversion/SCALE_NOW;
  double vel_rescale = vel_conversion*SCALE_NOW;
  float ds, dx, z, scale_now, z1, hubble;
  if (LIGHTCONE && SCALE_NOW != 1.0) Error("SCALE_NOW value inconsistent with LIGHTCONE\n");
  for (i=0; i<nelems; i++) {
      dx = 0.0f;
      for (j=0; j<3; j++) {
	  p[p_start+i].pos[j]   *= pos_rescale;
	  p[p_start+i].pos[j]   += 0.5*BOX_SIZE;
	  p[p_start+i].pos[j+3] *= vel_rescale;
	  ds = p[p_start+i].pos[j]-LIGHTCONE_ORIGIN[j];
	  dx += ds*ds;      
      }
      if (LIGHTCONE) { /* Remove hubble flow and convert from physical to comoving */
	  z = comoving_distance_h_to_redshift(sqrt(dx));
	  scale_now = scale_factor(z);
	  z1 = 1.0/scale_now;
	  hubble = 100.0*hubble_scaling(z1-1.0);
	  for (j=0; j<3; j++) {
	      p[p_start+i].pos[j+3] -= hubble * (p[p_start+i].pos[j]-LIGHTCONE_ORIGIN[j]);
	      p[p_start+i].pos[j+3] *= scale_now;
	  }
      }
  }
}

void load_particles_sdf(char *name, char *header, struct particle **p, int64_t *num_p)
{
    int i;
    SDF *sdfp;
    char filename[256];
    int num_blocks = NUM_BLOCKS;
    int block = 0;
    int64_t gnobj;
    int nobj;
    int64_t start;
    float particle_mass[2];
    char *mnames[1] = {"mass"};
    int nnames;
    int no_ident = 0;
    void *addrs[7];
    char *names[7] = {"x", "y", "z", "vx", "vy", "vz", "ident"};
    int strides[7];
    int nobjs[7];
    int64_t starts[7];
    int ret;
    double H0, R0, redshift, epsilon_scaled, epsilon;
    int units_rockstar = 0;

    sscanf(name, "%s %d", filename, &block);
    sdfp = SDFopen(header, filename);
    if (sdfp == NULL) singlError("[Error] SDFopen failed, %s", SDFerrstring);

    if (block == 0) fprintf(stderr, "libSDF Version: %s\n", libSDF_version);

    /* SDFget does not write third arg if key is not found */
    int fileversion = 0;
    SDFgetint(sdfp, "version", &fileversion);
    SDFgetint(sdfp, "version_2HOT", &fileversion);
    int units_2HOT = 0;
    SDFgetint(sdfp, "units_2HOT", &units_2HOT);
    if (fileversion == 1) {
	units_2HOT = 1;
	SDFgetdoubleOrDie(sdfp, "Omega0",  &Om);
	SDFgetdoubleOrDie(sdfp, "Lambda_prime",  &Ol);
	SDFgetdoubleOrDie(sdfp, "H0",  &H0);
	h0 = H0*10.0*(one_kpc/one_Gyr);
    } else if (fileversion == 2) {
	double Or, Of;
	units_2HOT = 1;
	SDFgetdouble(sdfp, "Omega0_m",  &Om);
	SDFgetdouble(sdfp, "Omega0_r",  &Or);
	SDFgetdouble(sdfp, "Omega0_lambda",  &Ol);
	SDFgetdouble(sdfp, "Omega0_fld",  &Of);
	Ol += Of;
	SDFgetdouble(sdfp, "h_100",  &h0);
	Om += Or;	/* put Or into Om, not strictly correct */
    } else {
	SDFgetint(sdfp, "units_rockstar", &units_rockstar);
	if (!units_rockstar) {
	    if (block == 0) fprintf(stderr, "[Warning] File version not specified, assuming rockstar units\n");
	    units_rockstar = 1;
	    SDFgetdouble(sdfp, "BOX_SIZE",  &BOX_SIZE);
	    SDFgetdouble(sdfp, "SCALE_NOW",  &SCALE_NOW);
	    SDFgetdouble(sdfp, "Om",  &Om);
	    SDFgetdouble(sdfp, "Omega0_m",  &Om);
	    SDFgetdouble(sdfp, "Ol",  &Ol);
	    SDFgetdouble(sdfp, "Omega0_lambda",  &Ol);
	    SDFgetdouble(sdfp, "h0",  &h0);
	    SDFgetdouble(sdfp, "H0",  &h0);
	}
    }
    if (fabs(Om + Ol - 1.0) > 1e-5) 
	singlError("[Error] Halo Finder Not Currently Configured to Run on Cosmologies with Curvature. (Omega_Matter = %f, Omega_Lambda = %f!)\n", Om, Ol)

    SDFgetint64(sdfp, "PERIODIC",  &PERIODIC);
    SDFgetint64(sdfp, "do_periodic",  &PERIODIC);

    if (SDFhasname("R0", sdfp)) {
	SDFgetdouble(sdfp, "R0",  &R0);
	BOX_SIZE = 2.0*R0;
    }

    if (SDFhasname("Rx", sdfp)) {
	SDFgetdouble(sdfp, "Rx",  &R0);
	BOX_SIZE = 2.0*R0;
    }

    if (SDFhasname("redshift", sdfp)) {
	SDFgetdouble(sdfp, "redshift", &redshift);
	SCALE_NOW = 1.0/(1.0+redshift);
    }

    if (SDFhasname("eps", sdfp)) { /* tree16 format */
	SDFgetdouble(sdfp, "eps", &epsilon);
	FORCE_RES = epsilon;
	FORCE_RES_PHYS_MAX = epsilon*SCALE_NOW;
    }
    if (SDFhasname("epsilon_scaled", sdfp)) {
	SDFgetdouble(sdfp, "epsilon_scaled", &epsilon_scaled);
	FORCE_RES = epsilon_scaled/SCALE_NOW;
	FORCE_RES_PHYS_MAX = epsilon_scaled;
    }

    SDFgetint64(sdfp, "light_cone", &LIGHTCONE);
    SDFgetdouble(sdfp, "light_cone_x0", &LIGHTCONE_ORIGIN[0]);
    SDFgetdouble(sdfp, "light_cone_y0", &LIGHTCONE_ORIGIN[1]);
    SDFgetdouble(sdfp, "light_cone_z0", &LIGHTCONE_ORIGIN[2]);

    if (units_2HOT) {
	BOX_SIZE *= h0/1000.0;
	FORCE_RES *= h0/1000.0;
	FORCE_RES_PHYS_MAX *= h0/1000.0;
	for (i=0; i<3; i++) {
	    LIGHTCONE_ORIGIN[i] *= h0/1000.0;
	    LIGHTCONE_ORIGIN[i] += BOX_SIZE/2.0;
	}
    }

    if (SDFgetint64(sdfp, "npart", &gnobj)) {
	gnobj = SDFnrecs(names[0], sdfp);
	if (block == 0) fprintf(stderr, "[Warning] %s does not have an \"npart\".\n", filename);
	if (block == 0) fprintf(stderr, "[Warning] Guessing npart=%ld from SDFnrecs(., %s)\n", 
		gnobj, names[0]);
    }
    TOTAL_PARTICLES = gnobj;

    if (SDFhasname(mnames[0], sdfp)) {
	nnames = 1;
	starts[0] = 0;
	nobjs[0] = 2;
	addrs[0] = &particle_mass[0];
	strides[0] = sizeof(float);
	ret = SDFseekrdvecsarr(sdfp, nnames, mnames, starts, nobjs, addrs, strides);
	if (ret != 0) Error("[Error] SDFread failed, %s\n", SDFerrstring);
	
	if (particle_mass[0] != particle_mass[1]) Error("[Error] unequal particle masses not supported\n");
	PARTICLE_MASS = particle_mass[0];
	if (units_2HOT) PARTICLE_MASS *= h0*1e10;
    } else if (SDFhasname("part_mass", sdfp)) {
	SDFgetdouble(sdfp, "part_mass", &PARTICLE_MASS);
    } else {
	SDFgetdouble(sdfp, "particle_mass", &PARTICLE_MASS);
	if (units_2HOT) PARTICLE_MASS *= h0*1e10;
    }

    if (block == 0) fprintf(stderr, "SDFread %s PARTICLE_MASS %g BOX_SIZE %f SCALE_NOW %f\nOm %f h0 %f FORCE_RES %f TOTAL_PARTICLES %ld\n", 
			    filename, PARTICLE_MASS, BOX_SIZE, SCALE_NOW, Om, h0, FORCE_RES, TOTAL_PARTICLES);
    if (PERIODIC) {
	double tmp = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3.0) / TOTAL_PARTICLES;
	/* universal constants differ, or would be 1e-5 */
	if (fabs(PARTICLE_MASS/tmp-1.0) > 1e-3) {
	    if (block == 0) { singlError("[Error] PARTICLE_MASS is not consistent with cosmological parameters " \
				    "%.4lg  %.4lg\n", PARTICLE_MASS, tmp);
	    } else exit(1);
	}
    }

    AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
    
    sdf_split(gnobj, num_blocks, block, &nobj, &start);

    *p = (struct particle *)check_realloc(*p, ((*num_p)+nobj)*sizeof(struct particle), "Allocating particles.");

    if (!SDFhasname("ident", sdfp)) no_ident = 1;
    if (no_ident) nnames = 6;
    else nnames = 7;
    for (i = 0; i < nnames; i++) {
	if (!SDFhasname(names[i], sdfp)) singlError("[Error] SDF file does not have %s\n", names[i]);
	starts[i] = start;
	nobjs[i] = nobj;
	strides[i] = sizeof(struct particle);
    }
    addrs[0] = &(*p)[0].pos[0];
    addrs[1] = &(*p)[0].pos[1];
    addrs[2] = &(*p)[0].pos[2];
    addrs[3] = &(*p)[0].pos[3];
    addrs[4] = &(*p)[0].pos[4];
    addrs[5] = &(*p)[0].pos[5];
    addrs[6] = &(*p)[0].id;
    ret = SDFseekrdvecsarr(sdfp, nnames, names, starts, nobjs, addrs, strides);
    if (ret != 0) Error("[Error] SDFread failed, %s\n", SDFerrstring);
    if (no_ident) {
	if (block == 0) fprintf(stderr, "[Warning] No ident in SDF file, numbering sequentially\n");
	sdf_renumber_particles(*p, *num_p, nobj, start);
    }
    sdf_check(*p, *num_p, nobj);
    if (units_2HOT) sdf_rescale_particles(*p, *num_p, nobj, h0/1000.0, one_kpc/one_Gyr);
    *num_p += nobj;

    SDFclose(sdfp);
}
