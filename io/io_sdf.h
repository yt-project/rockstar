#ifndef _IO_SDF_H_
#define _IO_SDF_H_
#include <stdint.h>
#include "../particle.h"

extern char libSDF_version[];
void load_particles_sdf(char *filename, char *header, struct particle **p, int64_t *num_p);

#endif /* _IO_SDF_H_ */
