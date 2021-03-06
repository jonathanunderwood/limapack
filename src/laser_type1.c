#include <stdlib.h>
#include <math.h>

#include "laser_type1.h"
#include "laser.h"
#include "memory.h"
#include "config_gsl_complex.h"
#include "au.h"

static laser_polzn_vector_t *
laser_type1_get_polzn_vector(const laser_t * self, const double t)
/* Returns the laser polarization vector. For type 1 lasers, this is
   independent of t. */
{
  laser_type1_t *l = (laser_type1_t *) self;
  return l->e;
} 

static double
laser_type1_get_envelope (const laser_t * self, const double t)
/* Returns the envelope function at time t. For a type 1 laser this is
   a pulse with a Gaussian profile, and is real valued. */
{
  laser_type1_t *l = (laser_type1_t *) self;
  double env;

  if (t < l->t0)
    {
      double a = (t - l->t0);
      double trise = l->trise;
      env = exp (-(a * a) / (2.0 * trise * trise));
    }
  else if (t > l->t0)
    {
      double a = (t - l->t0);
      double tfall = l->tfall;
      env = exp (-(a * a) / (2.0 * tfall * tfall));
    }
  else
    env = 1.0;

  /* envmin is the minimum value of the envelope function when it's
     normalized to 1 that we wish to consider. If env is less than
     this value we return a negative value signifying we wish to skip
     this laser at this time t. */
  if (env > l->envmin) 
    return l->E0 * env; 
  else return -1.0; 
}

static double
laser_type1_get_frequency (const laser_t * self, const double t)
/* Returns the laser frequency at time t. For type 1 lasers this is
   independent of t. */
{
  laser_type1_t *l = (laser_type1_t *) self;

  return l->freq;
}


static void 
laser_type1_dtor(laser_t *laser)
{
  laser_type1_t * l = (laser_type1_t *) laser;
  laser_polzn_vector_dtor(l->e);
  MEMORY_FREE(l);
}

laser_t *
laser_type1_cfg_parse_ctor (config_setting_t * element)
{
  double I, freq, t0, trise, tfall;
  laser_type1_t * laser;

  if (MEMORY_ALLOC(laser) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  /* Register dispatch functions. */
  laser_dispatch_register ((laser_t *) laser, 
			   &laser_type1_get_polzn_vector, 
			   &laser_type1_get_envelope, 
			   &laser_type1_get_frequency,
			   &laser_type1_dtor);

  if (!(config_setting_lookup_float (element, "t0", &t0) &&
	config_setting_lookup_float (element, "trise", &trise) &&
	config_setting_lookup_float (element, "tfall", &tfall) &&
	config_setting_lookup_float (element, "envmin", &(laser->envmin)) &&
	config_setting_lookup_float (element, "I", &I) &&
	config_setting_lookup_float (element, "freq", &freq) &&
	config_setting_lookup_float (element, "phi", &(laser->phi)) &&
	config_setting_lookup_float (element, "theta", &(laser->theta)) &&
	config_setting_lookup_float (element, "chi", &(laser->chi)) &&
	config_setting_lookup_gsl_complex (element, "ex", &(laser->ex)) &&
	config_setting_lookup_gsl_complex (element, "ey", &(laser->ey)) &&
	config_setting_lookup_gsl_complex (element, "ez", &(laser->ez))
	))
    {
      fprintf(stderr, "Incomplete laser type 1 description in file.\n"); 
      MEMORY_FREE (laser);
      return NULL;
    }

  /* Convert laser intensity in W/cm^2 to electric field in AU. */
  laser->E0 = E_TO_AU (LASER_I_TO_E(I));

  /* Convert frequency in THz to angular frequency in AU. */
  laser->freq = HZ_TO_AU (freq * 1.0e12) * 2.0 * M_PI;

  /* Convert laser temporal paraters from picoseconds to AU. */
  laser->t0 = PS_TO_AU (t0);
  laser->trise = PS_TO_AU (trise);
  laser->tfall = PS_TO_AU (tfall);

  /* Construct laser polarization vector in the lab frame. */
  laser->e = laser_polzn_vector_ctor_from_cart (laser->ex, laser->ey, laser->ez);
  laser->e->rotate(laser->e, laser->phi, laser->theta, laser->chi);

  return (laser_t *)laser;
}

