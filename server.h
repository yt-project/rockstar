#ifndef _SERVER_H_
#define _SERVER_H_

#define READER_TYPE 0
#define WRITER_TYPE 1
#define DATA_PARTICLES 0
#define DATA_HALOS 1

#include <stdint.h>

struct client_info {
  int type;
  int64_t cs;
  char *address, *serv_port;
  int port;
  float bounds[6];
  float bounds_prevsnap[6];
  int64_t num_halos;
  int64_t status;
  int64_t extra_workers;
};

#define PROJECTION_SIZE 10000

struct projection_request {
  int64_t dir, id;
  float bounds[6];
};

struct projection {
  int64_t dir, id;
  float bounds[6];
  int64_t data[PROJECTION_SIZE];
};


int server(void);
void check_num_writers(void);

#endif /* _SERVER_H_ */
