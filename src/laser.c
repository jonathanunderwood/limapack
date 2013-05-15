#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <libconfig.h>

#include "memory.h"
#include "laser.h"
#include "laser_polzn.h"
#include "laser_type1.h"

/* laser_t * */
/* laser_ctor (laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, const double t), */
/* 	    double (*get_envelope) (const laser_t *self, const double t), */
/* 	    double (*get_frequency) (const laser_t *self, const double t), */
/* 	    void (*dtor) (laser_t *self) */
/* 	    ) */
/* { */
/*   laser_t *l;  */

/*   if (MEMORY_ALLOC(l) < 0) */
/*     { */
/*       MEMORY_OOMERR; */
/*       return NULL; */
/*     } */

/*   laser_dispatch_register (l, get_polzn_vector, get_envelope, get_frequency, dtor);  */

/*   return l; */
/* } */

void
laser_dispatch_register (laser_t * laser,
			 laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, const double t),
			 double (*get_envelope) (const laser_t *self, const double t),
			 double (*get_frequency) (const laser_t *self, const double t),
			 void (*dtor) (laser_t *self)
			 )		
{
  laser->get_polzn_vector = get_polzn_vector;
  laser->get_envelope = get_envelope;
  laser->get_frequency = get_frequency;
  laser->dtor = dtor;
}

/* void  */
/* laser_dtor (laser_t *self) */
/* /\* This should actually never be called, but rather, the registered dispatch */
/*    function dtor should be. *\/ */
/* { */
/*   MEMORY_FREE(self); */
/* } */

laser_collection_t * 
laser_collection_ctor (const int nlasers)
{
  laser_collection_t * l;

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

void laser_collection_dtor (laser_collection_t * lasers)
{
  int i;

  for (i = 0; i < lasers->nlasers; i++)
    {
      /* We may not have completely filled the array of pointers, so some
	 elements could be still NULL - this would happen following a failure in
	 config parsing. */
      if (lasers->laser[i] != NULL)
	{
	  //((laser_t *)(lasers->laser[i]))->dtor((laser_t *)lasers->laser[i]);
	  lasers->laser[i]->dtor(lasers->laser[i]);
	}
    }

  MEMORY_FREE(lasers->laser);
  MEMORY_FREE(lasers);
}

laser_collection_t *
laser_collection_cfg_parse_ctor (const config_t * cfg)
{
  config_setting_t *setting;
  laser_collection_t *lasers = NULL;

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

      lasers = laser_collection_ctor(nlasers);
      if (lasers == NULL)
	return NULL;
      
      for(i = 0; i < nlasers; i++)
	{
	  config_setting_t *this_laser = config_setting_get_elem(setting, i);
	  laser_type_t type;
	  int ret;
	  
	  if (this_laser == NULL)
	    {
	      fprintf(stderr, "Failed to find a laser in the lasers section.\n");
	      laser_collection_dtor(lasers);
	      return NULL;
	    }

	  ret = config_setting_lookup_int(this_laser, "type", (long int *)(&type)) ;
	  
	  if (ret == CONFIG_FALSE)
	    {
	      fprintf(stderr, "Failed to get laser type.\n");
	      laser_collection_dtor(lasers);
	      return NULL;
	    }

	  switch (type) 
	    {
	    case LASER_TYPE1:
	      lasers->laser[i] = laser_type1_cfg_parse_ctor(this_laser);

	      if (lasers->laser[i] == NULL)
		{
		  laser_collection_dtor(lasers);
		  fprintf(stderr, "Error parsing type 1 laser configuration.\n");
		  return NULL;
		}
	      break;
		
	    default:
	      fprintf(stderr, "Unrecognized laser type: %d\n", type);
	      laser_collection_dtor(lasers);
	      return NULL;
	    }
	}

      return lasers;
    }
}

int 
laser_collection_all_negligible (laser_collection_t * lasers, const double t)
/* Returns 1 if all laser intensities are negligible at this time t, 0 otherwise. */
{
  int i;

  for (i = 0; i < lasers->nlasers; i++)
    {
      /* Our convention is that the get_envelope dispatch function
	 will return a value of less than zero if this laser intensity
	 is negligible. */
      if (lasers->laser[i]->get_envelope(lasers->laser[i], t) > 0.0)
	return 0;
    }

  return 1;

}
