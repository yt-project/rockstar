#ifndef _LITEHASH_H_
#define _LITEHASH_H_
#include <inttypes.h>

struct litebucket {
  void *key, *data;
  int64_t next;
};

struct litehash {
  uint64_t keywidth, hashwidth, elems, num_buckets, extra_buckets, hashnum,
    extra_buckets_alloced;
  struct litebucket *buckets;
};


struct litehash *new_litehash(uint64_t keywidth);
void *lh_keylist(struct litehash *lh);
void *lh_getval(struct litehash *lh, void *key);
void lh_setval(struct litehash *lh, void *key, void *data);
void free_litehash(struct litehash *lh);
void lh_prealloc(struct litehash *lh, int64_t size);

void lh_setval2(struct litehash *lh, void *key1, void *key2, void *data);
void *lh_getval2(struct litehash *lh, void *key1, void *key2);
void free_litehash2(struct litehash *lh);

#endif /* _LITEHASH_H_ */

