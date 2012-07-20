#include laser.h

void
laser_ctor (laser_t *self,
	    laser_polzn_vector * (*get_polzn_vector) (laser_t *self, const double t)
	    laser_polzn_tensor * (*get_polzn_tensor) (laser_t *self, const double t)
	    gsl_comlplex * (*get_envelope) (laser_t *self, const double t)
	    )
{
  self->get_polzn_vector = &get_polzn_vector;
  self->get_polzn_tensor = &get_polzn_tensor;
  self->get_envelope = &get_envelope;
}

void 
laser_xtor (laser_t *self)
{
  self->get_polzn_vector = NULL;
  self->get_polzn_tensor = NULL;
  self->get_envelope = NULL;
}
