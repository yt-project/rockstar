#ifndef CHECK_SYSCALLS_H
#define CHECK_SYSCALLS_H

void system_error(char *errmsg);
FILE *check_fopen(char *filename, char *mode);
FILE *check_popen(char *command, char *mode);
void *check_realloc(void *ptr, size_t size, char *reason);
size_t check_fread(void *ptr, size_t size, size_t nitems, FILE *stream);
size_t check_fwrite(void *ptr, size_t size, size_t nitems, FILE *stream);
void check_fseeko(FILE *stream, off_t offset, int whence);
char *check_fgets(char *ptr, size_t size, FILE *stream);
FILE *check_rw_socket(char *command, pid_t *pid);
void rw_socket_close(FILE *res, pid_t pid);

#define check_fprintf(file, ...) { if (fprintf(file, __VA_ARGS__) <= 0)	{  \
      fprintf(stderr, "[Error] Failed printf to fileno %d!\n", fileno(file)); \
      perror("[Error] Reason"); \
      exit(1); \
    }}

#endif /* CHECK_SYSCALLS_H */
