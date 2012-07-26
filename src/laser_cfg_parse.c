#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "laser.h"
#include "laser_cfg_parse.h"
#include "laser_type1.h"
#include "memory.h"

laser_container_t * 
laser_container_ctor (const int nlasers)
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

  for (i = 0; i < lasers->nlasers; i++)
    {
      /* We may not have completely filled the array of pointers, so some
	 elements could be still NULL - this would happen following a failure in
	 config parsing. */
      if (lasers->laser[i] != NULL)
	{
	  //lasers->laser[i] = (laser_t *)lasers->laser[i];
	  ((laser_t *)(lasers->laser[i]))->dtor((laser_t *)lasers->laser[i]);
	}
    }

  MEMORY_FREE(lasers->laser);
  MEMORY_FREE(lasers);
}

laser_container_t *
laser_container_cfg_parse_ctor (const config_t * cfg)
{
  config_setting_t *setting;
  laser_container_t *lasers = NULL;

  setting = config_lookup(cfg, "lasers");

  if (setting == NULL)
    {
      fprintf(stderr, "Failed to find lasers section in config.\n");
      return NULL;
    }
  else
    {
      int i;
      int nlasers = config_setting_length(setting);

      lasers = laser_container_ctor(nlasers);
      if (lasers == NULL)
	return NULL;
      
      for(i = 0; i < nlasers; i++)
	{
	  config_setting_t *this_laser = config_setting_get_elem(setting, i);
	  laser_type_t type = NONE;
	  int ret = config_setting_lookup_int(this_laser, "type", (int *)(&type)) ;
	  
	  if (ret == CONFIG_FALSE)
	    {
	      fprintf(stderr, "Failed to get laser type.\n");
	      laser_container_dtor(lasers);
	      return NULL;
	    }

	  switch (type) 
	    {
	    case TYPE1:
	      lasers->laser[i] = (laser_type1_t *) lasers->laser[i];
	      lasers->laser[i] = laser_type1_ctor();
	      if (lasers->laser[i] == NULL)
		{
		  laser_container_dtor(lasers);
		  return NULL;
		}

		ret = laser_type1_cfg_parse (this_laser, lasers->laser[i]);
		if (ret)
		  {
		    laser_container_dtor(lasers);
		    return NULL;
		  }
		
		break;
		
	    default:
	      fprintf(stderr, "Unrecognized laser type: %d\n", type);
	      laser_container_dtor(lasers);
	      return NULL;
	    }
	}

      return lasers;
    }
}
