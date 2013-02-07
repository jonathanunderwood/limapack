#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#include "memory.h"
#include "slurp.h"

int
slurp_file_to_buffer(const char *filename, char **buffer, const int max_file_size)
/* Slurp a file into a buffer and return it. max_file_size is the
   maximum file size on bytes. Pass a negative value to disable
   checking file is smaller than max_file_size. The return value is
   the size of the buffer. */
{
  struct stat file_status;
  FILE *fp = NULL;
  int buffer_size;
  
  if (stat(filename, &file_status) != 0)
    {
      perror("Could not stat or file does not exist:");
      return -1;
    }

  if ((max_file_size > 0) && (file_status.st_size > max_file_size))
    {
      fprintf(stderr, "Config file %s too large, size = %d kB\n.", 
	      filename, (int) (file_status.st_size / 1024));
      return -1;
    }
  
  buffer_size = file_status.st_size + 1; /* +1 for \0 at end */
  if (MEMORY_ALLOC_N((*buffer), buffer_size) < 0) 
    {
      MEMORY_OOMERR;
      return -1;
    }
    
  fp = fopen (filename, "r");
  if (fp == NULL)
    {
      fprintf(stderr, "Failed to open file %s. Exiting.", filename);
      perror("System error message:");
      MEMORY_FREE((*buffer));
      return -1;
    }

  if (!fread(*buffer, file_status.st_size, 1, fp)) 
    {
      fprintf(stderr, "Could not read file %s. Exiting.", filename);
      perror("System error message:");
      MEMORY_FREE (buffer);
      return -1;
    }
  fclose (fp);
  
  /* Add terminal null - fread doesn't. */
  (*buffer)[buffer_size - 1] = '\0';

  return buffer_size + 1;
}

/* For posterity, commented out below is an alternative strategy which
   doesn't rely on stat which is probably more portable. */
    
/* char * */
/* slurp_file (const char filename[]) */
/* { */
/*   char *buffer; */
/*   off_t size; */
/*   FILE *fnp = fopen (filename, "r"); */

/*   if (fnp == NULL) */
/*     return NULL; */

/*   fseek (fnp, 0, SEEK_END); */
/*   size = ftell (fnp); */

/*   /\* Add room for terminal NULL in case file doesn't end with a newline. *\/ */
/*   buffer = CHKMALLOC ((size + 1) * sizeof (char)); */

/*   rewind (fnp); */
/*   fread (buffer, size, 1, fnp); */
/*   fclose (fnp); */

/*   buffer[size] = '\0'; */

/*   return buffer; */
/* } */
