#include <stdlib.h>
#include <gsl/gsl_complex.h>

#include "memory.h"
#include "laser.h"
#include "laser_polzn.h"

laser_t *
laser_ctor (laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, const double t),
	    double (*get_envelope) (const laser_t *self, const double t),
	    double (*get_frequency) (const laser_t *self, const double t),
	    void (*dtor) (laser_t *self)
	    )
{
  laser_t *l; 

  if (MEMORY_ALLOC(l) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  laser_dispatch_register (l, get_polzn_vector, get_envelope, get_frequency, dtor); 

  return l;
}

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

void 
laser_dtor (laser_t *self)
/* This should actually never be called, but rather, the registered dispatch
   function dtor should be. */
{
  MEMORY_FREE(self);
}
