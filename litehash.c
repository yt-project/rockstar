#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include "check_syscalls.h"
#include "litehash.h"

#define MAX_LOAD_FACTOR 0.7

struct litehash *new_litehash(uint64_t keywidth) {
  struct litehash *lh = check_realloc(NULL, sizeof(struct litehash),
				       "Allocating litehash.");
  int64_t i;
  memset(lh, 0, sizeof(struct litehash));
  lh->keywidth = keywidth;
  lh->hashnum = rand() + 
    (((uint64_t)rand())<<(uint64_t)32) + (uint64_t)(rand());
  lh->hashwidth = 8;
  lh->num_buckets = (uint64_t)1 << lh->hashwidth;
  lh->buckets = check_realloc(NULL, sizeof(struct litebucket)*lh->num_buckets,
			      "Allocating hash buckets.");
  memset(lh->buckets, 0, sizeof(struct litebucket)*lh->num_buckets);
  for (i=0; i<lh->num_buckets; i++) lh->buckets[i].next = -1;
  if (!(lh->hashnum & 1)) lh->hashnum++;
  return lh;
}

void *lh_keylist(struct litehash *lh) {
  int64_t i, j=0;
  void *kl = check_realloc(NULL, lh->keywidth*lh->elems, "Allocating key list.");
  struct litebucket *lb = NULL;
  for (i=0; i<lh->num_buckets+lh->extra_buckets; i++) {
    lb = lh->buckets + i;
    if (!lb->key) continue;
    memcpy(kl + lh->keywidth*j, lb->key, lh->keywidth);
    j++;
  }
  return kl;
}


inline uint64_t _lh_hash_function(struct litehash *lh, void *key) {
  uint64_t intkey=0;
  if (lh->keywidth == sizeof(uint64_t)) intkey = *((uint64_t *)key);
  else if (lh->keywidth == sizeof(uint32_t)) intkey = *((uint32_t *)key);
  else memcpy(&intkey, key, lh->keywidth);
  return ((intkey*lh->hashnum)>>(64 - lh->hashwidth));
}

void _lh_add_more_buckets(struct litehash *lh, int64_t add_factor) {
  int64_t i;
  int64_t old_num_buckets = lh->num_buckets, old_extra_buckets = lh->extra_buckets;
  int64_t num_extra_buckets_filled = 0;
  struct litebucket *new_buckets, *lb, *newb;

  lh->num_buckets <<= add_factor;
  int64_t new_alloc_size = sizeof(struct litebucket)*(lh->num_buckets + lh->extra_buckets);

  new_buckets = check_realloc(NULL, new_alloc_size, "Allocating new buckets.");
  memset(new_buckets, 0, new_alloc_size);
  for (i=0; i<(lh->num_buckets+lh->extra_buckets); i++) new_buckets[i].next = -1;

  lh->hashwidth+=add_factor;
  for (i=0; i<old_num_buckets+old_extra_buckets; i++) {
    lb = lh->buckets + i;
    if (!lb->key) continue;
    newb = new_buckets + _lh_hash_function(lh, lb->key);
    if (newb->key) {
      new_buckets[lh->num_buckets + num_extra_buckets_filled] = *newb;
      newb->next = lh->num_buckets + num_extra_buckets_filled;
      num_extra_buckets_filled++;
    }
    newb->key = lb->key;
    newb->data = lb->data;
  }

  lh->extra_buckets_alloced = old_extra_buckets;
  lh->extra_buckets = num_extra_buckets_filled;
  free(lh->buckets);
  lh->buckets = new_buckets;
}


void lh_prealloc(struct litehash *lh, int64_t size) {
  if (size <= 0) return;
  int64_t numbits = ceil(log(size/MAX_LOAD_FACTOR)/log(2));
  if (numbits <= lh->hashwidth) return;
  _lh_add_more_buckets(lh, numbits-lh->hashwidth);
}

inline struct litebucket *_lh_getval(struct litehash *lh, void *key) {
  struct litebucket *lb = lh->buckets + _lh_hash_function(lh, key);
  while (lb->key) {
    if (lh->keywidth != sizeof(uint64_t)) {
      if (!memcmp(key, lb->key, lh->keywidth)) return lb;
    }
    else if (*((uint64_t *)key) == *((uint64_t *)lb->key)) return lb;
    if (lb->next < 0) return NULL;
    lb = lh->buckets + lb->next;
  }
  return NULL;
}

void *lh_getval(struct litehash *lh, void *key) {
  struct litebucket *lb = _lh_getval(lh, key);
  if (lb) return lb->data;
  return NULL;
}

void lh_setval(struct litehash *lh, void *key, void *data) {
  struct litebucket *lb;
  lb = _lh_getval(lh, key);
  if (!lb) {
    if (lh->elems>=lh->num_buckets*MAX_LOAD_FACTOR) _lh_add_more_buckets(lh,1);
    lb = lh->buckets + _lh_hash_function(lh, key);
    if (lb->key) {
      if (lh->extra_buckets >= lh->extra_buckets_alloced) {
	lh->extra_buckets_alloced += 1000;
	lh->buckets = check_realloc(lh->buckets, sizeof(struct litebucket)*
				    (lh->extra_buckets_alloced+lh->num_buckets),
				    "Allocating extra buckets.");
	lb = lh->buckets + _lh_hash_function(lh, key);
      }
      lh->buckets[lh->num_buckets+lh->extra_buckets] = *lb;
      lb->next = lh->num_buckets+lh->extra_buckets;
      lh->extra_buckets++;
    }
    lh->elems++;
    lb->key = key;
  }
  lb->data = data;
}

void free_litehash(struct litehash *lh) {
  if (!lh) return;
  if (lh->buckets) {
    memset(lh->buckets, 0, sizeof(struct litebucket)*
	   (lh->extra_buckets+lh->num_buckets));
    free(lh->buckets);
  }
  memset(lh, 0, sizeof(struct litehash));
  free(lh);
}

void free_litehash2(struct litehash *lh) {
  int64_t i;
  if (!lh) return;
  for (i=0; i<lh->num_buckets + lh->extra_buckets; i++)
    if (lh->buckets[i].key) free_litehash(lh->buckets[i].data);
  free_litehash(lh);
}


void lh_setval2(struct litehash *lh, void *key1, void *key2, void *data) {
  struct litehash *lh2 = lh_getval(lh, key1);
  if (!lh2) {
    lh2 = new_litehash(lh->keywidth);
    lh_setval(lh, key1, lh2);
  }
  lh_setval(lh2, key2, data);
}

void *lh_getval2(struct litehash *lh, void *key1, void *key2) {
  struct litehash *lh2 = lh_getval(lh, key1);
  if (!lh2) return NULL;
  return lh_getval(lh2, key2);
}
