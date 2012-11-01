#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "molecule.h"
#include "laser.h"
#include "odesys.h"
#include "memory.h"

int
main(int argc, char *argv[])
{
  config_t *cfg = NULL;
  molecule_t *molecule = NULL;
  laser_collection_t *lasers = NULL;
  odesys_t *odesys = NULL;
  int ret = 0;

  if (argc < 2)
    {
      fprintf (stderr, "Useage: %s <parameter file> <output file>.\n",
               argv[0]);
      exit (1);
    }

  /* Parse configurtion file into a libconfig config_t object. */ 
  if (MEMORY_ALLOC(cfg) < 0)
    {
      fprintf(stderr, "Failed to allocate memory for config object.\n");
      exit (1);
    }
  config_init (cfg);

  if (config_read_file (cfg, argv[1]) == CONFIG_FALSE)
    {
      switch (config_error_type(cfg))
	{
	case CONFIG_ERR_FILE_IO:
	  fprintf(stderr, "Failed to open file %s.\n", argv[1]);
	  ret = 1; 
	  goto cleanup_exit;
	  break;
	case CONFIG_ERR_PARSE:
	  fprintf(stderr, "Error in configuration file %s at line %d\n",
		  config_error_file(cfg), config_error_line(cfg));
	  fprintf(stderr, "Error was: \"%s\"\n", 
	      config_error_text(cfg));
	  ret = 1; 
	  goto cleanup_exit;
	  break;
	default:
	  /* do nothing. */
	  break;
	}
    }
  /* Parse the config_t object. */
  molecule = molecule_cfg_parse_ctor (cfg);
  if (molecule == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse molecule information from configuration file.\n");
      ret = 1; 
      goto cleanup_exit;
    }

  lasers = laser_collection_cfg_parse_ctor (cfg);
  if (lasers == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse laser information from configuration file.\n");
      ret = 1; 
      goto cleanup_exit;
    }
  
  odesys = odesys_cfg_parse_ctor (cfg);
  if (odesys == NULL)
    {
      fprintf(stderr, 
	      "Failed to initialize odesys.\n");
      ret = 1; 
      goto cleanup_exit;
    }

  config_destroy (cfg);
  free (cfg);
  cfg = NULL;

  odesys_init (odesys, molecule, lasers);

  odesys_tdse_propagate_simple (odesys);

 cleanup_exit:
  if (odesys != NULL)
    odesys_dtor(odesys);
  if (lasers != NULL)
    laser_collection_dtor(lasers);
  if (molecule != NULL)
    molecule_dtor (molecule);
  if (cfg != NULL)
    {
      config_destroy (cfg);
      MEMORY_FREE(cfg);
    }
  return ret;
}
