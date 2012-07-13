#ifndef _IO_BGC2_H_
#define _IO_BGC2_H_

#include <stdint.h>

extern char **bgc2_snapnames;
extern int64_t num_bgc2_snaps;

void output_bgc2(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds);
void calc_bgc2_parents(int64_t snap);

#endif /* _IO_BGC2_H_ */
