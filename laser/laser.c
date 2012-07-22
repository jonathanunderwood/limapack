#include "laser.h"

laser_t *
laser_ctor (laser_polzn_vector * (*get_polzn_vector) (const laser_t *self, const double t)
	    gsl_comlplex * (*get_envelope) (const laser_t *self, const double t)
	    double (*get_frequency) (const laser_t *self, const double t)
	    )
{
  laser_t *l; 

  if (MEMORY_ALLOC(l) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  laser_callback_register (l, get_polzn_vector, get_envelope, get_frequency); 

  return l;
}

void
laser_callback_register (laser_t * laser,
			 laser_polzn_vector * (*get_polzn_vector) (const laser_t *self, const double t)
			 gsl_comlplex * (*get_envelope) (const laser_t *self, const double t)
			 double (*get_frequency) (const laser_t *self, const double t)
			 )		
{
  laser->get_polzn_vector = get_polzn_vector;
  laser->get_envelope = get_envelope;
  laser->get_frequency = get_frequency;
}

void 
laser_xtor (laser_t *self)
{
  self->get_polzn_vector = NULL;
  self->get_polzn_tensor = NULL;
  self->get_envelope = NULL;
}
