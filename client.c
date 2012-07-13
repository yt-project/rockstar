#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <sys/select.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>
#include "check_syscalls.h"
#include "particle.h"
#include "rockstar.h"
#include "groupies.h"
#include "inet/rsocket.h"
#include "io/meta_io.h"
#include "config_vars.h"
#include "server.h"
#include "client.h"
#include "merger.h"
#include "distance.h"
#include "bounds.h"
#include "fun_times.h"
#include "universe_time.h"

#define CLIENT_DEBUG 0

struct recipient *recipients = NULL;
int64_t num_recipients = 0;
int64_t *part_ids = NULL;
int64_t *part_id_buffer = NULL;
int64_t num_buffer_part_ids = 0;

struct projection *prj = NULL;
struct projection_request *prq = NULL;
int64_t num_proj = 0;
int64_t in_error_state = 0;
int64_t RECIPIENT_BUFFER=100000;

void network_io_err(int64_t s) {
  if (in_error_state) return;
  in_error_state = 1;
  send_to_socket(s, "err!", 4);
}

void network_error_cleanup() {
  p = check_realloc(p, 0, "Freeing particle memory");
  num_p = 0;
  halos = check_realloc(halos, 0, "Freeing halos");
  num_halos = 0;
  rockstar_cleanup();
  clear_merger_tree();
  in_error_state = 0;
}

void reset_projection_count(void) {
  prj = check_realloc(prj, sizeof(struct projection)*num_proj,
		      "Allocating projections.");
  prq = check_realloc(prq, sizeof(struct projection_request)*num_proj,
		      "Allocating projection requests.");
}

struct recipient *add_recipient(int halo_recipient, int64_t c) {
  struct recipient *r;
  int64_t struct_size;
  num_recipients++;
  recipients = check_realloc(recipients, (sizeof(struct recipient)
					  * num_recipients), 
			     "Allocating particle data recipients.");
  r = recipients+num_recipients-1;
  memset(r, 0, sizeof(struct recipient));
  struct_size = (halo_recipient) ? sizeof(struct halo):sizeof(struct particle);
  r->buffer = check_realloc(NULL, struct_size*RECIPIENT_BUFFER,
			    "Allocating recipient transmit buffer.");

  recv_from_socket(c, r->bounds, sizeof(float)*6);
  r->address = recv_msg_nolength(c, r->address);
  r->port = recv_msg_nolength(c, r->port);
  return r;
}

void clear_recipients(void) {
  int64_t i;
  for (i=0; i<num_recipients; i++) {
    free(recipients[i].buffer);
    free(recipients[i].port);
    free(recipients[i].address);
  }
  free(recipients);
  recipients = NULL;
  num_recipients = 0;
}

void calc_particle_bounds(float *bounds) {
  int64_t i,j;
  for (j=0; j<6; j++) bounds[j]=0;
  if (!num_p) return;
  memcpy(bounds, p[0].pos, sizeof(float)*3);
  memcpy(bounds+3, p[0].pos, sizeof(float)*3);
  for (i=1; i<num_p; i++) {
    for (j=0; j<3; j++) {
      if (bounds[j] > p[i].pos[j]) 
	bounds[j] = p[i].pos[j];
      if (bounds[j+3] < p[i].pos[j]) 
	bounds[j+3] = p[i].pos[j];
    }
  }
}

void trim_particles(float *bounds) {
  int64_t i, j;
  if (!TRIM_OVERLAP) return;
  for (i=0; i<3; i++) {
    bounds[i] += TRIM_OVERLAP;
    bounds[i+3] -= TRIM_OVERLAP;
  }

  if (ROUND_AFTER_TRIM)
    for (i=0; i<6; i++) 
      bounds[i]=((int64_t)(bounds[i]/ROUND_AFTER_TRIM + 0.5))*ROUND_AFTER_TRIM;

  for (i=0; i<num_p; i++)
    for (j=0; j<3; j++)
      if (bounds[j] > p[i].pos[j] || bounds[j+3] <= p[i].pos[j]) {
	num_p--;
	p[i] = p[num_p];
	i--;
	break;
      }

  p = check_realloc(p, sizeof(struct particle)*num_p, "Removing overlap.");
}

void clear_particle_rbuffer(struct recipient *r) {
  if (!r->buffered) return;
  send_to_socket(r->cs, "part", 4);
  send_msg(r->cs, r->buffer, sizeof(struct particle)*r->buffered);
  r->buffered = 0;
}

void clear_halo_rbuffer(struct recipient *r) {
  int64_t i, pids=0;
  struct halo *bh = r->buffer;
  if (!r->buffered) return;
  send_to_socket(r->cs, "halo", 4);
  send_msg(r->cs, r->buffer, sizeof(struct halo)*r->buffered);

  for (i=0; i<r->buffered; i++) pids+=bh[i].num_p;
  if (pids>num_buffer_part_ids) {
    part_id_buffer = check_realloc(part_id_buffer, sizeof(int64_t)*pids,
				   "Allocating particle ID buffer.");
    num_buffer_part_ids = pids;
  }
  pids = 0;
  for (i=0; i<r->buffered; i++) {
    memcpy(part_id_buffer + pids, part_ids + bh[i].p_start,
	   sizeof(int64_t)*bh[i].num_p);
    pids += bh[i].num_p;
  }
  send_to_socket(r->cs, "pids", 4);
  send_msg(r->cs, part_id_buffer, sizeof(int64_t)*pids);
  r->buffered = 0;
}


void add_particle_to_buffer(struct recipient *r, struct particle *p1) {
  struct particle *buffer = r->buffer;
  if (r->buffered == RECIPIENT_BUFFER) clear_particle_rbuffer(r); 
  buffer[r->buffered] = *p1;
  r->buffered++;
}

void add_halo_to_buffer(struct recipient *r, struct halo *h1) {
  struct halo *buffer = r->buffer;
  if (r->buffered == RECIPIENT_BUFFER) clear_halo_rbuffer(r); 
  buffer[r->buffered] = *h1;
  r->buffered++;
}

void check_particle_bounds(struct particle *p1, struct recipient *r) {
  struct particle pt = *p1;
  if (_check_bounds(p1->pos, pt.pos, r->bounds)) add_particle_to_buffer(r, &pt);
}

void check_halo_bounds(struct halo *h1, struct recipient *r) {
  struct halo ht = *h1;
  if (_check_bounds(h1->pos, ht.pos, r->bounds)) add_halo_to_buffer(r, &ht);
}

int64_t check_projection_bounds(struct particle *p1, struct projection *pr) {
  return (_check_bounds_raw(p1->pos, pr->bounds));
}

void send_config(int64_t c) {
  double data;
#define snd(x) { data = x; send_to_socket_noconfirm(c, &data, sizeof(double)); }
  snd(PARTICLE_MASS);
  snd(AVG_PARTICLE_SPACING);
  snd(SCALE_NOW);
  snd(BOX_SIZE);
  snd(Ol);
  snd(Om);
  snd(h0);
#undef snd
}

void recv_config(int64_t c) {
  int64_t i;
#define rcv(x) { recv_from_socket(c, &x, sizeof(double)); }
  rcv(PARTICLE_MASS);
  rcv(AVG_PARTICLE_SPACING);
  rcv(SCALE_NOW);
  rcv(BOX_SIZE);
  rcv(Ol);
  rcv(Om);
  rcv(h0);
#undef rcv
  for (i=0; i<3; i++)
    if (LIGHTCONE_ORIGIN[i] || LIGHTCONE_ALT_ORIGIN[i]) break;
  if (i==3) {
    for (i=0; i<3; i++) LIGHTCONE_ORIGIN[i] = LIGHTCONE_ALT_ORIGIN[i] = BOX_SIZE/2.0;
  }
}

void send_particles(int64_t c) {
  int64_t i,j;
  for (j=0; j<num_recipients; j++) 
    recipients[j].cs = 
      connect_to_addr(recipients[j].address, recipients[j].port);
  for (i=num_p-1; i>=0; i--) {
    for (j=0; j<num_recipients; j++) check_particle_bounds(p+i, recipients+j);
    if (!(i%PARTICLE_REALLOC_NUM)) 
      p = check_realloc(p,sizeof(struct particle)*i,"Freeing particle memory.");
  }
  num_p = 0;
  p = check_realloc(p,0,"Freeing particle memory.");
  for (j=0; j<num_recipients; j++) {
    clear_particle_rbuffer(recipients+j);
    send_to_socket(recipients[j].cs, "done", 4);
    close_rsocket(recipients[j].cs);
  }
  clear_recipients();
  send_to_socket(c, "done", 4);
}

void send_halos(char *c_address, char *c_port, int64_t snap, int64_t chunk) {
  int64_t i,j, c;
  struct binary_output_header bheader;

  load_binary_halos(snap, chunk, &bheader, &halos, &part_ids);

  for (j=0; j<num_recipients; j++)
    recipients[j].cs=connect_to_addr(recipients[j].address, recipients[j].port);

  for (i=0; i<bheader.num_halos; i++)
    for (j=0; j<num_recipients; j++) check_halo_bounds(halos+i, recipients+j);

  for (j=0; j<num_recipients; j++) {
    clear_halo_rbuffer(recipients+j);
    send_to_socket(recipients[j].cs, "cnfg", 4);
    send_config(recipients[j].cs);
    send_to_socket(recipients[j].cs, "done", 4);
    close_rsocket(recipients[j].cs);
  }
  
  c = connect_to_addr(c_address, c_port);
  send_to_socket(c, "rdne", 4);
  exit(0);
}

void close_connection(int64_t cs, int64_t *cslist, int64_t *num_cs) {
  int64_t i;
  for (i=0; i<*num_cs; i++) {
    if (cslist[i] == cs) {
      close_rsocket(cslist[i]);
      *num_cs = (*num_cs)-1;
      cslist[i] = cslist[*num_cs];
      break;
    }
  }
}

void receive_stuff(int64_t s, int64_t c, int64_t timestep) {
  int64_t i, j, num_senders = 0, done = 0, length, new_p_start, *senders = NULL;
  char cmd[5] = {0};
  struct binary_output_header *bheader;
  struct halo **halos_recv, *th;
  int64_t max_conn, **pids_recv;

  while (!in_error_state && (num_senders || !done)) {
    clear_rsocket_tags();
    if (!done) {
      tag_rsocket(s);
      tag_rsocket(c);
    }
    for (i=0; i<num_senders; i++) tag_rsocket(senders[i]);
    max_conn = select_rsocket(RSOCKET_READ, 0);
    for (i=0; i<max_conn; i++) {
      if (!check_rsocket_tag(i)) continue;
      if (i==s) {
	num_senders++;
	senders = check_realloc(senders, sizeof(int64_t)*num_senders,
				"Allocating particle sender FDs.");
	senders[num_senders-1] = accept_connection(s,NULL,NULL);
      }
      else if (i==c) {
	recv_from_socket(c, cmd, 4);
	if (!strcmp(cmd, "done")) { done = 1; }
	else if (!strcmp(cmd, "err!")) {
	  in_error_state = 1;
	  for (j=0; j<num_senders; j++) close_rsocket(senders[j]);
	  num_senders = 0;
	  break;
	}
	else { fprintf(stderr, "[Error] Server protocol error rs (%s)!\n", cmd); exit(1); }
      }
      else {
	if (recv_from_socket(i, cmd, 4)<=0) {
	  network_io_err(c);
	  for (j=0; j<num_senders; j++) close_rsocket(senders[j]);
	  num_senders = 0;
	  break;
	}
	if (!strcmp(cmd, "part")) {
	  length = num_p*sizeof(struct particle);
	  p = recv_msg(i, p, &length, length);
	  assert(!(length%(sizeof(struct particle))));
	  num_p = length / sizeof(struct particle);
	}

	else if (!strcmp(cmd, "halo")) {
	  assert(timestep > 0);
	  bheader = (timestep > 1) ? &head2 : &head1;
	  pids_recv = (timestep > 1) ? &part2 : &part1;
	  halos_recv = (timestep > 1) ? &halos2 : &halos1;

	  length = bheader->num_halos*sizeof(struct halo);
	  *halos_recv = recv_msg(i, *halos_recv, &length, length);
	  assert(!(length%sizeof(struct halo)));
	  length /= sizeof(struct halo);

	  //Redo particle pointers
	  new_p_start = bheader->num_particles;
	  for (j=bheader->num_halos; j<length; j++) {
	    th = (*halos_recv) + j;
	    th->p_start = new_p_start;
	    new_p_start += th->num_p;
	  }
	  bheader->num_halos = length;

	  recv_from_socket(i, cmd, 4);
	  assert(!strcmp(cmd, "pids"));
	  length = bheader->num_particles*sizeof(int64_t);
	  *pids_recv = recv_msg(i, *pids_recv, &length, length);
	  assert(!(length%sizeof(int64_t)));
	  bheader->num_particles = length / sizeof(int64_t);
	}

	else if (!strcmp(cmd, "cnfg")) {
	  recv_config(i);
	}

	else if (!strcmp(cmd, "done"))
	  close_connection(i, senders, &num_senders);

	else if (!strcmp(cmd, "rdne")) {
	  close_connection(i, senders, &num_senders);
	  send_to_socket_noconfirm(c, "done", 4);
	}

	else { fprintf(stderr, "[Error] Client protocol error rs (%s)!\n", cmd); exit(1); }
      }
    }
  }
  free(senders);
}

void do_projections(void) {
  int64_t i, j, idx, dir;
  assert(BOX_SIZE > 0);
  for (i=0; i<num_proj; i++) {
    prj[i].id = prq[i].id;
    dir = prj[i].dir = prq[i].dir;
    memcpy(prj[i].bounds, prq[i].bounds, sizeof(float)*6);
    for (j=0; j<PROJECTION_SIZE; j++) prj[i].data[j] = 0;
    for (j=0; j<num_p; j++) {
      if (check_projection_bounds(p+j, prj+i)) {
	idx = (double)PROJECTION_SIZE*p[j].pos[dir]/(double)BOX_SIZE;
	if (idx >= PROJECTION_SIZE) idx = PROJECTION_SIZE-1;
	prj[i].data[idx]++;
      }
    }
  }
}

void accept_workloads(char *c_address, char *c_port, int64_t snap, int64_t chunk) {
  char *address = NULL, *port = NULL;
  int64_t s = connect_to_addr(c_address, c_port);
  int64_t m = s;
  int64_t id = -1, new_bounds = 1;
  char cmd[5] = {0};
  float zero_bounds[6];
  struct fof *fofs = NULL;
  struct workunit_info w;
  rockstar_cleanup();
  p = check_realloc(p, 0, "Freeing particles.\n");
  num_p = 0;

  if (s < 0) exit(1);
  memset(zero_bounds, 0, sizeof(float)*6);

  send_to_socket(m, &chunk, sizeof(int64_t));
  while (1) {
    if (recv_from_socket(m, cmd, 4)<=0) {
      fprintf(stderr, "[Warning] Failed to receive instruction from server for chunk %"PRId64" (shutting down)!\n", chunk);
      send_to_socket(s, "clos", 4);
      exit(1);
    }
    if (!strcmp(cmd, "wrku")) {
      recv_from_socket(m, &w, sizeof(struct workunit_info));
      assert(w.num_particles >= 0 && w.num_fofs >= 0 && !w.num_halos);
      num_p = w.num_particles;
      fofs = recv_and_alloc(m, fofs, sizeof(struct fof)*w.num_fofs);
      p = recv_and_alloc(m, p, sizeof(struct particle)*w.num_particles);
      if (CLIENT_DEBUG) fprintf(stderr, "Received %"PRId64" particles and %"PRId64" fofs from id %"PRId64" (Worker %"PRId64")\n", w.num_particles, w.num_fofs, id-NUM_READERS, chunk);
      if (new_bounds && TEMPORAL_HALO_FINDING) {
	new_bounds = 0;
	if (!memcmp(w.bounds, zero_bounds, sizeof(float)*6))
	  load_previous_halos(snap, w.chunk, NULL);
	else load_previous_halos(snap, w.chunk, w.bounds);
      }
      do_workunit(&w, fofs);
      if (CLIENT_DEBUG) fprintf(stderr, "Analyzed %"PRId64" particles and %"PRId64" fofs from id %"PRId64", and found %"PRId64" halos. (Worker %"PRId64")\n", w.num_particles, w.num_fofs, id-NUM_READERS, num_halos, chunk);
      w.num_halos = num_halos;
      send_to_socket(m, "wrkd", 4);
      send_to_socket(m, &w, sizeof(struct workunit_info));
      send_to_socket(m, fofs, sizeof(struct fof)*w.num_fofs);
      send_to_socket(m, halos, sizeof(struct halo)*w.num_halos);
      send_to_socket(m, extra_info, sizeof(struct extra_halo_info)*w.num_halos);
      send_to_socket(m, p, sizeof(struct particle)*w.num_particles);
    }
    else if (!strcmp(cmd, "nmwk")) {
      if (m != s) close_rsocket(m);
      p = check_realloc(p, 0, "Freeing particles.");
      send_to_socket(s, "nmwk", 4);
      send_to_socket(s, &id, sizeof(int64_t));
      m = s;
    }
    else if (!strcmp(cmd, "work")) {
      assert(m==s);
      recv_from_socket(s, &id, sizeof(int64_t));
      address = recv_msg_nolength(s, address);
      port = recv_msg_nolength(s, port);
      m = connect_to_addr(address, port);
      if (m<0) {
	send_to_socket(s, "clos", 4);
	exit(1);
      }
      send_to_socket(m, &chunk, sizeof(int64_t));
      new_bounds = 1;
    }
    else if (!strcmp(cmd, "fini")) {
      exit(0);
    }
    else {
      fprintf(stderr, "[Error] Error in client protocol aw (%s)!\n", cmd);
      send_to_socket_noconfirm(s, "clos", 4);
      exit(1);
    }
  }
}

int send_workunit(int64_t sock, struct workunit_info *w,
		  struct fof **fofs, struct particle **particles,
		  int64_t *no_more_work) {
  if (!(*no_more_work)) find_unfinished_workunit(w, fofs, particles);
  if ((*no_more_work) || (!w->num_fofs)) {
    *no_more_work = 1;
    send_to_socket(sock, "nmwk", 4);
    return 0;
  }
  send_to_socket(sock, "wrku", 4);
  send_to_socket(sock, w, sizeof(struct workunit_info));
  send_to_socket(sock, *fofs, sizeof(struct fof)*w->num_fofs);
  send_to_socket(sock, *particles, sizeof(struct particle)*w->num_particles);
  if ((*particles >= p) && (*particles < p+num_p)) *particles = NULL;
  return 1;
}

void distribute_workloads(int64_t c, int64_t s, int64_t snap, int64_t chunk, float *bounds) {
  int64_t i, j, num_workers = 0, no_more_work = 0, workdone = 0,
    done = 0, id_offset, worker_chunk, hcnt, id;
  int64_t *workers = NULL, max_cs = 0, child = -1, new_w, child_has_connected=0;
  char *address = NULL, *port = NULL;
  char cmd[5] = {0};
  struct workunit_info w;
  struct particle *parts = NULL;
  struct fof *fofs = NULL;
  struct halo *rhalos = NULL;
  struct extra_halo_info *ehi = NULL;

  if (bounds) memcpy(w.bounds, bounds, sizeof(float)*6);
  else memset(w.bounds, 0, sizeof(float)*6);
  w.chunk = chunk;

  while ((num_workers || !done) && !(in_error_state)) {

    if (work_finished() && !workdone && no_more_work && child_has_connected) {
      workdone = 1;
      rockstar_cleanup();
      hcnt = count_halos_to_print(bounds);
      send_to_socket_noconfirm(c, "hcnt", 4);
      send_to_socket_noconfirm(c, &hcnt, sizeof(int64_t));
      if (CLIENT_DEBUG) fprintf(stderr, "Analysis of %"PRId64" halos complete (chunk %"PRId64")\n", hcnt, chunk);
    }

    clear_rsocket_tags();
    if (!done || child > -1) tag_rsocket(c);
    if (!done) tag_rsocket(s);
    for (i=0; i<num_workers; i++) tag_rsocket(workers[i]);

    max_cs = select_rsocket(RSOCKET_READ, 0);
    for (i=0; i<max_cs; i++) {
      if (!check_rsocket_tag(i)) continue;
      if (i==s) {
	new_w = accept_connection(s,NULL,NULL);
	recv_from_socket(new_w, &worker_chunk, sizeof(int64_t));
	if (worker_chunk == chunk) {
	  assert(child < 0);
	  child = new_w;
	  child_has_connected = 1;
	}
	if (send_workunit(new_w, &w, &fofs, &parts, &no_more_work) ||
	    (child == new_w)) {
	  workers = check_realloc(workers, sizeof(int64_t)*(num_workers+1),
				  "Allocating rockstar analysis FDs.");
	  workers[num_workers] = new_w;
	  num_workers++;
	  if (CLIENT_DEBUG) fprintf(stderr, "Got new worker (chunk %"PRId64"; wchunk %"PRId64")\n", chunk, worker_chunk);
	} else {
	  close_rsocket(new_w);
	}
      }
      else if (i==c) {
	recv_from_socket(c, cmd, 4);
	if (!strcmp(cmd, "fini")) {
	  assert(child >= 0);
	  send_to_socket(child, "fini", 4);
	  close_connection(child, workers, &num_workers);
	  child = -1;
	  if (CLIENT_DEBUG) fprintf(stderr, "Told child to finish (chunk %"PRId64")\n", chunk);
	}
	else if (!strcmp(cmd, "work")) {
	  assert(child >= 0);
	  recv_from_socket(c, &id, sizeof(int64_t));
	  address = recv_msg_nolength(c, address);
	  port = recv_msg_nolength(c, port);
	  send_to_socket(child, "work", 4);
	  send_to_socket(child, &id, sizeof(int64_t));
	  send_msg(child, address, strlen(address)+1);
	  send_msg(child, port, strlen(port)+1);
	  if (CLIENT_DEBUG) fprintf(stderr, "Child (%"PRId64") connecting to id %"PRId64"\n", chunk, id);
	}
	else if (!strcmp(cmd, "outp")) {
	  assert(no_more_work && !done);
	  recv_from_socket(c, &id_offset, sizeof(int64_t));
	  output_and_free_halos(id_offset, snap, chunk, bounds);
	  p = check_realloc(p, 0, "Freeing particle memory");
	  num_p = 0;
	  done = 1;
	  if (CLIENT_DEBUG) fprintf(stderr, "Finished (chunk %"PRId64")\n", chunk);
	}
	else if (!strcmp(cmd, "quit")) {
	  exit(1);
	}
	else if (!strcmp(cmd, "err!")) {
	  in_error_state = 1;
	  for (j=0; j<num_workers; j++) close_rsocket(workers[j]);
	  num_workers = 0;
	  break;
	}
	else { fprintf(stderr, "[Error] Server protocol error dw (%s)!\n", cmd); exit(1); }
      }
      else {
	if (recv_from_socket(i, cmd, 4) <= 0) {
	  network_io_err(c);
	  for (j=0; j<num_workers; j++) close_rsocket(workers[j]);
	  num_workers = 0;
	  break;
	}
	if (!strcmp(cmd, "wrkd")) {
	  recv_from_socket(i, &w, sizeof(struct workunit_info));
	  assert(w.num_particles >= 0 && w.num_particles <= num_p);
	  fofs = recv_and_alloc(i, fofs, sizeof(struct fof)*w.num_fofs);
	  rhalos = recv_and_alloc(i, rhalos, sizeof(struct halo)*w.num_halos);
	  ehi = recv_and_alloc(i, ehi, sizeof(struct extra_halo_info)*w.num_halos);
	  parts = recv_and_alloc(i, parts, sizeof(struct particle)*w.num_particles);
	  integrate_finished_workunit(&w, fofs, rhalos, ehi, parts);
	  if (!send_workunit(i, &w, &fofs, &parts, &no_more_work) && 
	      (child != i))
	    close_connection(i, workers, &num_workers);
	}
	else if (!strcmp(cmd, "clos")) {
	  close_connection(i, workers, &num_workers);
	}
	else if (!strcmp(cmd, "nmwk")) {
	  assert(i == child);
	  recv_from_socket(i, &id, sizeof(int64_t));
	  send_to_socket_noconfirm(c, "nmwk", 4);
	  send_to_socket_noconfirm(c, &id, sizeof(int64_t));
	}

	else { fprintf(stderr, "[Error] Client protocol error dw (%s)!\n", cmd); exit(1); }
      }
    }
  }
  check_realloc(address, 0, "Freeing address memory.");
  check_realloc(port, 0, "Freeing port memory.");
  check_realloc(workers, 0, "Freeing worker memory.");
  check_realloc(fofs, 0, "Freeing FOF memory.");
  check_realloc(parts, 0, "Freeing particle memory.");
  check_realloc(rhalos, 0, "Freeing halo memory.");
  check_realloc(ehi, 0, "Freeing extra halo info memory.");
}


void client(void) {
  int64_t snap, block, chunk, type=-1, i, n, timestep, id_offset;
  struct binary_output_header bheader;
  uint64_t magic;
  char buffer[1024];
  char cmd[5] = {0};
  char *hostname = NULL;
  char port[10] = {0};
  float bounds[6], box_size;
  //struct recipient *r;
  int64_t c, s = -1, portnum, readers;
  int stat_loc;
  int64_t num_nodes;

  clear_merger_tree();

  if (FORK_READERS_FROM_WRITERS) {
    num_nodes = (FORK_PROCESSORS_PER_MACHINE) ? 
      (NUM_WRITERS/FORK_PROCESSORS_PER_MACHINE) : NUM_WRITERS;
    readers = NUM_READERS/num_nodes + ((NUM_READERS%num_nodes) ? 1 : 0);
    for (i=0; i<readers; i++) {
      n = fork();
      if (n==0) { type = READER_TYPE; break; } 
      if (n<0) system_error("Couldn't fork reader process!");
    }
    if (i==readers) { type = WRITER_TYPE; }
  }

  if (FORK_PROCESSORS_PER_MACHINE && (type != READER_TYPE)) {
    for (i=1; i<FORK_PROCESSORS_PER_MACHINE; i++) {
      n = fork();
      if (n==0) break;
      if (n<0) system_error("Couldn't fork process!");
    }
  }

  c = connect_to_addr(PARALLEL_IO_SERVER_ADDRESS, PARALLEL_IO_SERVER_PORT);
  recv_from_socket(c, &magic, sizeof(uint64_t));
  if (magic != ROCKSTAR_MAGIC) {
    fprintf(stderr, "[Error] Received invalid client responses.  Check network connectivity.\n");
    exit(1);
  }
  send_to_socket_noconfirm(c, &magic, sizeof(uint64_t));
  if (type < 0) send_to_socket_noconfirm(c, "rdwr", 4);
  if (type == READER_TYPE) send_to_socket_noconfirm(c, "read", 4);
  if (type == WRITER_TYPE) send_to_socket_noconfirm(c, "writ", 4);

  while (1) {
    recv_from_socket(c, cmd, 4);
    if (!strcmp(cmd, "err!")) network_error_cleanup();

    if (!strcmp(cmd, "read")) {
      type = READER_TYPE;
      hostname = recv_msg_nolength(c, hostname);
      send_msg(c, port, sizeof(char));
    }

    else if (!strcmp(cmd, "writ")) {
      type = WRITER_TYPE;
      hostname = recv_msg_nolength(c, hostname);
      for (i=0; i<5000 && s<0; i+=29) {
	portnum = PARALLEL_IO_WRITER_PORT+i;
	snprintf(port, 10, "%"PRId64, portnum);
	s = listen_at_addr(hostname, port);
      }
      if (i>=5000) {
	fprintf(stderr, "[Error] Couldn't start particle data server at %s:%d-%d!\n",
		hostname, (int)PARALLEL_IO_WRITER_PORT, (int)portnum);
	exit(1);
      }
      send_msg(c, port, strlen(port)+1);
      if (LIGHTCONE && strlen(LIGHTCONE_ALT_SNAPS)) {
	memcpy(LIGHTCONE_ORIGIN, LIGHTCONE_ALT_ORIGIN, sizeof(double)*3);
      }
      for (i=0; i<3; i++)
	if (LIGHTCONE_ORIGIN[i] || LIGHTCONE_ALT_ORIGIN[i]) break;
    }

    else if (!strcmp(cmd, "snap"))
      recv_from_socket(c, &snap, sizeof(int64_t));

    else if (!strcmp(cmd, "chnk"))
      recv_from_socket(c, &chunk, sizeof(int64_t));

    else if (!strcmp(cmd, "rdbk")) {
      assert(type == READER_TYPE);
      recv_from_socket(c, &block, sizeof(int64_t));
      if (LIGHTCONE && strlen(LIGHTCONE_ALT_SNAPS) && block >= (NUM_BLOCKS/2)) {
	if (LIGHTCONE == 1)
	  read_input_names(LIGHTCONE_ALT_SNAPS, &snapnames, &NUM_SNAPS);
	LIGHTCONE = 2;
	get_input_filename(buffer, 1024, snap, block-(NUM_BLOCKS/2));
      }
      else {
	if (LIGHTCONE == 2) {
	  LIGHTCONE = 1;
	  read_input_names(SNAPSHOT_NAMES, &snapnames, &NUM_SNAPS);
	}	  
	get_input_filename(buffer, 1024, snap, block);
      }
      read_particles(buffer);
    }
      
    else if (!strcmp(cmd, "cnf?")) {
      assert(type == READER_TYPE);
      calc_particle_bounds(bounds);
      if (TRIM_OVERLAP) trim_particles(bounds);
      send_to_socket(c, "bxsz", 4);
      box_size = BOX_SIZE;
      send_to_socket(c, &box_size, sizeof(float));
      send_to_socket(c, bounds, sizeof(float)*6);
      send_to_socket(c, "cnfg", 4);
      send_config(c);
    }

    else if (!strcmp(cmd, "cnfg")) {
      assert(type == WRITER_TYPE);
      recv_config(c);
      if (LIGHTCONE) init_cosmology();
      init_time_table();
    }

    else if (!strcmp(cmd, "rdy?")) {
      send_to_socket(c, "rdy!", 4);
    }

    else if (!strcmp(cmd, "proj")) {
      assert(type == READER_TYPE);
      recv_from_socket(c, &num_proj, sizeof(int64_t));
      reset_projection_count();
      recv_from_socket(c, prq, sizeof(struct projection_request)*num_proj);
      if (CLIENT_DEBUG) fprintf(stderr, "Doing client projections for block %"PRId64"\n", block);
      do_projections();
      if (CLIENT_DEBUG) fprintf(stderr, "Done with client projections for block %"PRId64"\n", block);
      send_to_socket(c, "cprj", 4);
      send_to_socket(c, &num_proj, sizeof(int64_t));
      for (i=0; i<num_proj; i++)
	send_to_socket(c, prj+i, sizeof(struct projection));
    }

    else if (!strcmp(cmd, "bnds")) {
      assert(type == WRITER_TYPE);
      recv_from_socket(c, bounds, sizeof(float)*6);
    }

    else if (!strcmp(cmd, "rcpt")) {
      assert(type == READER_TYPE);
      add_recipient(PARTICLE_RECIPIENT, c);
    }

    else if (!strcmp(cmd, "xfrp")) {
      if (type == READER_TYPE) send_particles(c);
      if (type == WRITER_TYPE) receive_stuff(s,c,0);
      if (in_error_state) network_error_cleanup();
    }

    else if (!strcmp(cmd, "rock")) {
      assert(type == WRITER_TYPE);
      recv_from_socket(c, &chunk, sizeof(int64_t));
      rockstar(bounds, 1);
      if (CLIENT_DEBUG) fprintf(stderr, "Found %"PRId64" fofs in chunk %"PRId64"\n", num_all_fofs, chunk);
      n = fork();
      if (!n) accept_workloads(hostname, port, snap, chunk);
      if (n<0) system_error("Couldn't fork halo analysis process!");
      distribute_workloads(c, s, snap, chunk, bounds);
      waitpid(n, &stat_loc, 0);
      if (in_error_state) network_error_cleanup();
      send_to_socket(c, "done", 4);
    }

    else if (!strcmp(cmd, "outp")) {
      assert(type == WRITER_TYPE);
      id_offset = 0; //Only used for dumping particles
      output_and_free_halos(id_offset, snap, chunk, bounds);
      p = check_realloc(p, 0, "Freeing particle memory");
      num_p = 0;
      send_to_socket(c, "done", 4);
    }

    else if (!strcmp(cmd, "gbds")) {
      load_binary_header(snap, chunk, &bheader);
      send_to_socket(c, &(bheader.bounds), sizeof(float)*6);
    }

    else if (!strcmp(cmd, "rcph")) {
      assert(type == WRITER_TYPE);
      add_recipient(HALO_RECIPIENT, c);
    }

    else if (!strcmp(cmd, "xfrh")) {
      assert(type == WRITER_TYPE);
      recv_from_socket(c, &timestep, sizeof(int64_t));
      n = fork();
      if (!n) send_halos(hostname, port, snap, chunk);
      if (n<0) system_error("Couldn't fork merger tree process!");
      receive_stuff(s,c,timestep);
      waitpid(n, &stat_loc, 0);
      if (in_error_state) {
	network_error_cleanup();
	continue;
      }
      if (timestep == 1) init_descendants();
      else if (timestep == 2) connect_particle_ids_to_halo_ids();
      clear_recipients();
    }

    else if (!strcmp(cmd, "merg")) {
      assert(type == WRITER_TYPE);
      calculate_descendants();
    }

    else if (!strcmp(cmd, "outc")) {
      assert(type == WRITER_TYPE);
      output_merger_catalog(snap, chunk, halos1, head1.num_halos);
      send_to_socket(c, "done", 4);
      clear_merger_tree();
    }

    else if (!strcmp(cmd, "delb")) {
      assert(type == WRITER_TYPE);
      if (DELETE_BINARY_OUTPUT_AFTER_FINISHED)
	delete_binary(snap, chunk);
    }
    
    else if (!strcmp(cmd, "quit")) {
      close_rsocket(c);
      return;
    }
  }
}
