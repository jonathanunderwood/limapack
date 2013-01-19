#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include <mpi.h>

#include "molecule.h"
#include "laser.h"
#include "odesys.h"
#include "memory.h"
#include "slurp.h"

#define MAX_CFG_FILE_SIZE (50*1024) // 50 kB

static int
parse_config_from_buffer(const char *buffer, molecule_t *molecule, 
			 laser_collection_t *lasers, odesys_t *odesys)
{
  config_t cfg;

  molecule = NULL;
  lasers = NULL;
  odesys = NULL;

  /* Parse configurtion file into a libconfig config_t object. */ 
  config_init (&cfg);

  if (config_read_string (&cfg, buffer) == CONFIG_FALSE)
    {
      switch (config_error_type(&cfg))
	{
	case CONFIG_ERR_PARSE:
	  fprintf(stderr, "Error in configuration at line %d\n",
		  config_error_line(&cfg));
	  fprintf(stderr, "Error was: \"%s\"\n", 
	      config_error_text(&cfg));
	  return -1; 
	default:
	  fprintf(stderr, "Unknown error in parsing configuration.\n");
	  return -1;
	}
    }

  /* Parse the config_t object. */
  molecule = molecule_cfg_parse_ctor (&cfg);
  if (molecule == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse molecule information from configuration file.\n");
      config_destroy (&cfg);
      return -1; 
    }

  lasers = laser_collection_cfg_parse_ctor (&cfg);
  if (lasers == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse laser information from configuration file.\n");
      molecule->dtor(molecule);
      config_destroy (&cfg);
      return -1; 
    }
  
  odesys = odesys_cfg_parse_ctor (&cfg);
  if (odesys == NULL)
    {
      fprintf(stderr, 
	      "Failed to initialize odesys.\n");
      molecule->dtor (molecule);
      laser_collection_dtor (lasers);
      config_destroy (&cfg);
      return -1; 
    }

  config_destroy (&cfg);

  return 0;
}

int
main(int argc, char *argv[])
{
  int ret = 0;
  char host[1024];
  int host_length, rank;
  molecule_t *molecule = NULL;
  laser_collection_t *lasers = NULL;
  odesys_t *odesys = NULL;
  int cfg_file_buff_size;
  char * cfg_file_buff;

  if (argc < 2)
    {
      fprintf (stderr, "Useage: %s <parameter file> <output file>.\n",
               argv[0]);
      exit (1);
    }

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(host, &host_length);

  fprintf(stdout, "<%d> Host: %s\n", rank, host);

  if (rank == 0) /* We are the parent process. */
    {
      /* molecule_t *molecule = NULL; */
      /* laser_collection_t *lasers = NULL; */
      /* odesys_t *odesys = NULL; */
      /* char *cfg_file_buff; */
      char *outfile;
      //      int cfg_file_buff_size;
      int i, nprocs;

      /* Read config file into a buffer for distribution to all processes
	 - this saves having to serialize. */
      cfg_file_buff_size = slurp_file_to_buffer(argv[1], &cfg_file_buff, MAX_CFG_FILE_SIZE);
      
      if (cfg_file_buff_size < 0)
	{
	  fprintf(stderr, "Failed to read in config file %s. Exiting.\n",
		  argv[1]);
	  goto exit;
	}

      /* Find number of child processes running. */
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      
      /* Distribute the parameter file to child processes. */
      for (i = 1; i < nprocs; i++)
	{
	  fprintf (stdout, "<%d> Sending parameters to process %d.\n", rank, i);
	  MPI_Send (&cfg_file_buff_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	  MPI_Send (cfg_file_buff, cfg_file_buff_size, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	}
      
      /* Parse config file - even though we're the parent process we
	 need to access some information. */
      if (parse_config_from_buffer(cfg_file_buff, molecule, lasers, odesys))
	goto exit;

      MEMORY_FREE (cfg_file_buff);

      odesys_init(odesys, molecule, lasers);
      
      odesys_tdse_propagate_mpi_master(odesys);
      
      outfile = argv[2];
      odesys_expval_fwrite(odesys, outfile);
      

    }
  else /* We are a slave child process. */
    {

      /* Receive config file buffer. */
      MPI_Recv (&cfg_file_buff_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (MEMORY_ALLOC_N(cfg_file_buff, cfg_file_buff_size) < 0)
	{
	  MEMORY_OOMERR;
	  goto exit;
	}

      MPI_Recv (cfg_file_buff, cfg_file_buff_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, 
		MPI_STATUS_IGNORE);
      
      if (parse_config_from_buffer(cfg_file_buff, molecule, lasers, odesys))
	goto exit;

      MEMORY_FREE (cfg_file_buff);

      odesys_init(odesys, molecule, lasers);
      
      odesys_tdse_propagate_mpi_slave(odesys);


    }


 exit:
  if (odesys != NULL)
    odesys_dtor(odesys);
  if (lasers != NULL)
    laser_collection_dtor(lasers);
  if (molecule != NULL)
    molecule->dtor (molecule);
  return ret;
}
