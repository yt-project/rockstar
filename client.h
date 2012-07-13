#ifndef _CLIENT_H_
#define _CLIENT_H_

#define PARTICLE_RECIPIENT 0
#define HALO_RECIPIENT 1
#define PARTICLE_REALLOC_NUM 100000


struct recipient {
  int64_t cs;
  char *address;
  char *port;
  float bounds[6];
  void *buffer;
  int64_t buffered;
};

void client(void);
void recv_config(int64_t c);
void send_config(int64_t c);

#endif /* _CLIENT_H_ */
