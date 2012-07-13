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
	$(CC) $(CFLAGS) load_bgc2.c ../check_syscalls.c ../io/io_util.c ../litehash.c -o load_bgc2  $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) load_full_particles.c ../check_syscalls.c  ../litehash.c ../io/stringparse.c -o load_full_particles  $(EXTRA_FLAGS)

clean:
	rm -f *~ 

