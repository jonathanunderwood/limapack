#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <getopt.h>

#include "odesys.h"

#define MAX_CONFIG_FILE_SIZE (64 * 1024) /* 64 kB */

static void
print_usage (const char *exec_name)
{
  printf ("Useage: %s [--nthreads=N] PARAM_FILE OUTPUT_FILE\n\n", exec_name); 
  printf ("\tN\tNumber of threads to employ (should be greater than 0)\n");
}

int
main(int argc, char *argv[])
{
  odesys_t *odesys = NULL;
  struct stat file_status;
  int nthreads = 1;
  char *paramfile, *outfile;

  /* Parse command line. */
  while (1)
    {
      int c;

      static struct option long_options[] =
	{
	  {"nthreads", required_argument, 0, 'n'},
	  {0, 0, 0, 0}
	};

      /* `getopt_long' stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "n:",
		       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
	break;

      switch (c)
	{
	case 'n':
	  sscanf (optarg, "%d", &nthreads);
	  break;
	case '?':
	  /* `getopt_long' already printed an error message. */
	  print_usage (argv[0]);
	  exit (1);
	  break;
	default:
	  abort ();
	}
    }

  if (nthreads < 0)
    {
      fprintf (stderr, "Error: nthreads should be >=0. Exiting.\n");
      exit (1);
    }

  if (argc - optind != 2)
    {
      print_usage (argv[0]);
      exit (1);
    }
  
  /* At this point, argv[optind] will be the first non-option
     argument. */
  paramfile = argv[optind];
  outfile = argv[optind + 1];

  if (stat(outfile, &file_status) == 0)
    {
      fprintf (stderr, "Output filename already exists. Exiting.\n");
      exit (1);
    }

  odesys = odesys_parse_cfg_from_file_ctor (paramfile, MAX_CONFIG_FILE_SIZE);
  if (odesys == NULL)
    {
      fprintf(stderr, "Failed to initialize. Exiting.\n");
      return -1;
    }

  if (nthreads == 1)
    odesys_tdse_propagate_simple (odesys);
  else
    odesys_tdse_propagate_threaded(odesys, nthreads);

  odesys_expval_fwrite(odesys, outfile);

  odesys_dtor(odesys);

  return 0;
}

#undef MAX_CONFIG_FILE_SIZE
