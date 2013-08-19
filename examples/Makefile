CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -lm -fno-math-errno -Wall
OFLAGS = -O3 -std=c99
DEBUGFLAGS = -g -O3 -std=c99
PROFFLAGS = -g -pg -O2 -std=c99
CC = gcc

all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROFFLAGS)"

reg:
	$(CC) $(CFLAGS) bgc2_to_ascii_particles.c load_bgc2.c ../check_syscalls.c ../io/io_util.c -o bgc2_to_ascii_particles  $(EXTRA_FLAGS)
	$(CC) -DTEST_LOADFP $(CFLAGS) load_full_particles.c ../check_syscalls.c   ../io/stringparse.c ../config_vars.c -o load_full_particles  $(EXTRA_FLAGS)
	$(CC) -DCALC_POTENTIALS $(CFLAGS) calc_potentials.c load_full_particles.c ../check_syscalls.c  ../io/stringparse.c ../config_vars.c ../potential.c -o calc_potentials  $(EXTRA_FLAGS)

clean:
	rm -f *~ 

