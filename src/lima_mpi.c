#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#include <mpi.h>

#include "odesys.h"
#include "slurp.h"
#include "memory.h"

#define MAX_CFG_FILE_SIZE (64*1024) // 64 kB
#define MAX_HOSTNAME_SIZE 1024

int
main(int argc, char *argv[])
{
  int ret = 0;
  char host[MAX_HOSTNAME_SIZE];
  int host_length, rank;
  odesys_t *odesys = NULL;
  int cfg_file_buff_size;
  char *cfg_file_buff;
  struct stat file_status;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(host, &host_length);

  fprintf(stdout, "<%d::%s Initiliasing.>\n", rank, host);

  if (argc < 2)
    {
      if (rank == 0)
	fprintf (stderr, "Useage: %s <parameter file> <output file>.\n",
		 argv[0]);
      MPI_Finalize();
      exit (-1);
    }

  /* Read config file if we're the master MPI process and distribute
     to other processes - this means that the config file does not
     need to be available on slave nodes. */
  if (rank == 0) /* Master process. */
    {
      if (stat(argv[2], &file_status) == 0)
	{
	  fprintf (stderr, "Output filename already exists. Exiting.\n");
	  MPI_Abort (MPI_COMM_WORLD, -1);
	  exit (-1);
	}

      /* Master process reads in file to a buffer and then broadcasts it. */
      cfg_file_buff_size = slurp_file_to_buffer(argv[1], &cfg_file_buff, 
						MAX_CFG_FILE_SIZE);
      
      if (cfg_file_buff_size < 0)
	{
	  fprintf(stderr, "Failed to read config file %s. Exiting.\n",
		  argv[1]);
	  MPI_Abort (MPI_COMM_WORLD, -1);
	  exit (-1);
	}
    }

  /* Distribute the config file buffer to all MPI slave processes. */
  MPI_Bcast(&cfg_file_buff_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank !=0)
    {
      if (MEMORY_ALLOC_N(cfg_file_buff, cfg_file_buff_size) < 0)
	{
	  MEMORY_OOMERR;
	  fprintf (stderr, "Failed to allocate storage for cfg_file_buff on %s. Slave exiting.\n",
		   host);
	  MPI_Abort (MPI_COMM_WORLD, -1);
	}
    }

  MPI_Bcast(cfg_file_buff, cfg_file_buff_size, MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Parse config file. */
  odesys = odesys_parse_cfg_from_buffer_ctor(cfg_file_buff);
  if (odesys == NULL)
    {
      fprintf (stderr, "Failed to parse config file %s.\n", argv[1]);
      MEMORY_FREE (cfg_file_buff);
      MPI_Abort (MPI_COMM_WORLD, -1);
    }

  MEMORY_FREE (cfg_file_buff);

  if (rank == 0) /* We are the parent process. */
    {
      char *outfile;
      
      if (odesys_tdse_propagate_mpi_master(odesys) < 0)
	{
	  fprintf(stderr, 
		 "odesys_tdse_propagate_mpi_master did not return cleanly.\n");
	  ret = -1;
	}
      else
	{
	  outfile = argv[2];
	  odesys_expval_fwrite(odesys, outfile);
	  ret = 0;
	}
      fprintf(stdout, "<%d::%s> Master exiting.\n", rank, host);
    }
  else /* We are a slave child process. */
    {
      if (odesys_tdse_propagate_mpi_slave(odesys) < 0)
	{
	  fprintf(stderr, 
		  "odesys_tdse_propagate_mpi_slave failed to complete cleanly.\n");
	  ret = -1;
	}
      else
	{
	  fprintf(stdout, "<%d::%s> Slave exiting.\n", rank, host);
	  ret = 0;
	}
    }

  odesys_dtor (odesys);
  
  MPI_Finalize ();

  exit (ret);
}

#undef MAX_CFG_FILE_SIZE
#undef MAX_HOSTNAME_SIZE
