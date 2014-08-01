CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall -fno-math-errno -fPIC $(SDF_INCLUDE)
LDFLAGS=-shared
OFLAGS = -lm -g -O3 -std=c99
DEBUGFLAGS = -lm -g -O0 -std=c99 -Dinline=
PROFFLAGS = -lm -g -pg -O2 -std=c99
CC = gcc
CFILES = rockstar.c check_syscalls.c fof.c groupies.c subhalo_metric.c potential.c nfw.c jacobi.c fun_times.c interleaving.c universe_time.c distance.c config_vars.c config.c bounds.c inthash.c io/read_config.c client.c server.c merger.c version.c inet/socket.c inet/rsocket.c inet/address.c io/meta_io.c io/io_internal.c io/io_ascii.c io/stringparse.c io/io_gadget.c io/io_generic.c io/io_art.c io/io_tipsy.c io/io_bgc2.c io/io_util.c
DIST_FLAGS =

HAS_SDF = $(shell pkg-config --exists libSDF-1.0 && echo "1" || echo "0")
ifeq ($(HAS_SDF), 1)
	SDF_INCLUDE = `pkg-config libSDF-1.0 --cflags`
	SDF_LIB = `pkg-config libSDF-1.0 --libs`
	CFLAGS += $(SDF_INCLUDE) -DENABLE_SDF
	OFLAGS += $(SDF_LIB)
	DEBUGFLAGS += $(SDF_LIB)
	PROFFLAGS += $(SDF_LIB)
	CFILES += io/io_sdf.c
# local SDF install:
# git clone https://bitbucket.org/JohnSalmon/SDF ; (cd SDF ; make)
else ifneq ("$(wildcard SDF/libSDF.a)","") 
	SDF_INCLUDE = -I./SDF
	SDF_LIB = -L./SDF -lSDF
	CFLAGS += $(SDF_INCLUDE) -DENABLE_SDF
	OFLAGS += $(SDF_LIB)
	DEBUGFLAGS += $(SDF_LIB)
	PROFFLAGS += $(SDF_LIB)
	CFILES += io/io_sdf.c
endif


all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)" EXT="-g"

prof:
	@make reg EXTRA_FLAGS="$(PROFFLAGS)" EXT="-p"

.REMAKE:

dist: .REMAKE
	cd ../ ; perl -ne 'print "$$1\n" if (/VERSION\s*\"([^\"]+)/)' Rockstar/version.h > Rockstar/VERSION; tar -czvf rockstar.tar.gz Rockstar/Makefile Rockstar/*.[ch] Rockstar/examples/Makefile Rockstar/[^sto]*/*.[ch] Rockstar/quickstart.cfg Rockstar/parallel.cfg Rockstar/scripts/*.pbs Rockstar/scripts/*.cfg Rockstar/scripts/*.pl Rockstar/SOURCE_LAYOUT Rockstar/README Rockstar/LICENSE Rockstar/VERSION Rockstar/ACKNOWLEDGMENTS Rockstar/CHANGELOG; mv rockstar.tar.gz Rockstar

versiondist:
	$(MAKE) dist DIST_FLAGS="$(DIST_FLAGS)"
	rm -rf dist
	mkdir dist
	cd dist; tar xzf ../rockstar.tar.gz ; perl -ne '/\#define.*VERSION\D*([\d\.]+)/ && print $$1' Rockstar/version.h > NUMBER ; mv Rockstar Rockstar-`cat NUMBER`; tar czf rockstar-`cat NUMBER`.tar.gz Rockstar-`cat NUMBER`

reg: version.c
	$(CC) $(CFLAGS) main.c $(CFILES) -o rockstar$(EXT)  $(EXTRA_FLAGS)

lib: version.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CFILES) -o librockstar.so $(OFLAGS)

bgc2:
	$(CC) $(CFLAGS) io/extra_bgc2.c mswutil/msw_bgc2.c $(CFILES) -o mswutil/msw_bgc2 $(OFLAGS)
	$(CC) $(CFLAGS) mswutil/msw_find_parents.c io/stringparse.c check_syscalls.c version.c -o mswutil/msw_find_parents -lm -O3 -std=c99 $(SDF_LIB)
	$(CC) $(CFLAGS) io/extra_bgc2.c util/redo_bgc2.c $(CFILES) -o util/finish_bgc2  $(OFLAGS)
	$(CC) $(CFLAGS) io/extra_bgc2.c util/bgc2_to_ascii.c $(CFILES) -o util/bgc2_to_ascii  $(OFLAGS)

parents:
	$(CC) $(CFLAGS) util/find_parents.c io/stringparse.c check_syscalls.c  -o util/find_parents $(OFLAGS)

substats:
	$(CC) $(CFLAGS) util/subhalo_stats.c $(CFILES) -o util/subhalo_stats  $(OFLAGS)

.PHONY: version.proto
version.proto:
	@echo char Rockstar_version\[\] = "\"`git describe --tags  --dirty` \"__DATE__\" \"__TIME__;" > version.proto

version.c: version.proto
	@cmp -s $< $@ || cp -p $< $@

clean:
	rm -f *~ io/*~ inet/*~ util/*~ rockstar util/redo_bgc2 util/subhalo_stats

