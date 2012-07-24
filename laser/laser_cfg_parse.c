#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "laser_type1.h"

laser_container_t * laser_container_ctor (const int nlasers)
{
  laser_container_t * l;

  if (MEMORY_ALLOC(l) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  if (MEMORY_ALLOC_N(l->laser, nlasers) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(l);
      return NULL;
    }

  l->nlasers = nlasers;

  return l;
}

void laser_container_dtor (laser_container_t * lasers)
{
  int i;

  for (i = 0; i<lasers->nlasers; i++)
    {
      /* We may not have completely filled the array of pointers, so some
	 elements could be still NULL. */
      if (lasers->laser[i] != NULL)
	MEMORY_FREE(lasers->laser[i]);
    }

  MEMORY_FREE(lasers->laser);
  MEMORY_FREE(lasers);
}


int
laser_cfg_parse (const config_t * cfg, laser_container_t *lasers)
{
  config_setting_t *s;

  s = config_lookup(&cfg, "lasers");

  if (setting == NULL)
    {
      fprintf(stderr, "Failed to find lasers section in config.\n");
      return -1;
    }
  else
    {
      int i;
      int nlasers = config_setting_length(setting);

      lasers = laser_container_ctor(nlasers);
      if (lasers == NULL)
	{
	  MEMORY_OOMERR;
	  return -1;
	}
      
      for(i = 0; i < nlasers; i++)
	{
	  config_setting_t *this_laser = config_setting_get_elem(setting, i);
	  laser_type_t type = NONE;
	  int ret = config_setting_lookup_int(this_laser, "type", &type) ;
	  
	  if (ret == CONFIG_FALSE)
	    {
	      fprintf(stderr, "Failed to get laser type at line %d:\n",
		      config_error_line (this_laser));
	      fprintf(stderr, "%s\n", config_error_text (this_laser));
	      laser_container_dtor(lasers);
	      return -1;
	    }

	  switch (type) 
	    {
	    case TYPE1:
	      if (MEMORY_ALLOC((laser_type1 *) (lasers->laser[i]) < 0))
		{
		  MEMORY_OOMERR;
		  laser_container_dtor(lasers);
		  return -1;
		}

		ret = laser_type1_parse (this_laser, lasers->laser[i]);
		if (ret)
		  {
		    laser_container_dtor(lasers);
		    return -1;
		  }
		
		break;
		
	    default:
	      fprintf(stderr, "Unrecognized laser type: %d\n", type);
	      laser_container_dtor(lasers);
	      return -1;
	    }
	}

      return 0;
    }
}
