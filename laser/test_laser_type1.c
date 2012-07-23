#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "laser_type1.h"

int
main ()
{
  config_t cfg;
  config_setting_t *setting;

  config_init (&cfg);

  if(! config_read_file(&cfg, "test_laser_type1.cfg"))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EXIT_FAILURE);
  }

  setting = config_lookup(&cfg, "lasers");

  if (setting == NULL)
    {
      fprintf(stderr, "Failed to find lasers section in file.\n");
      config_destroy(&cfg);
      exit (EXIT_FAILURE);
    }
  else
    {
      int nlasers = config_setting_length(setting);
      int i;
      laser_type1_t *laser_arr;

      if (MEMORY_ALLOC_N(laser_arr, nlasers) < 0) 
	{
	  MEMORY_OOMERR;
	  config_destroy(&cfg);
	  exit(1);
	}
      
      for(i = 0; i < nlasers; ++i)
	{
	  config_setting_t *this_laser = config_setting_get_elem(setting, i);
	  
	  int ret = laser_type1_parse (this_laser, &(laser_arr[i]));
	  
	  if (ret)
	    {
	      fprintf(stderr, "Failed to parse laser %d\n");
	      config_destroy(&cfg);
	      MEMORY_FREE (laser_arr);
	      exit (EXIT_FAILURE);
	    }
	}
    } 
  
  config_destroy(&cfg);

  MEMORY_FREE (laser_arr);
  exit (EXIT_SUCCESS);

}
