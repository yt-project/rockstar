#ifndef _STRINGPARSE_H_
#define _STRINGPARSE_H_

enum parsetype {
  PARSE_FLOAT32 = 0,
  PARSE_INT32 = 1,
  PARSE_FLOAT64 = 2,
  PARSE_INT64 = 3,
  PARSE_STRING = 4,
  PARSE_SKIP = 5,
};

#define SHORT_PARSETYPE enum short_parsetype {\
  F = 0, \
  D = 1, \
  F64 = 2, \
  LF = 2, \
  LD = 3, \
  D64 = 3, \
  S = 4, \
  K = 5 \
};



int stringparse(char *buffer, void **data, enum parsetype *types, int max_n);

#endif
