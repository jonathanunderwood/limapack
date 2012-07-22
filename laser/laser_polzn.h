#ifndef __LASER_POLZN_H__
#define __LASER_POLZN_H__ 1

#include <gsl/gsl_complex.h>

/* Polarization vector struct and functions. */

typedef struct _laser_polzn_vector
{
  gsl_complex (*get) (const struct _laser_polzn_vector * e, const int p);
  void (*set) (struct _laser_polzn_vector * e, const int p, gsl_complex val);
  void (*rotate) (struct _laser_polzn_vector * e,
		  const double phi, const double theta, const double chi);
  gsl_complex e[3];
} laser_polzn_vector_t;

laser_polzn_vector_t *
laser_polzn_vector_ctor_from_cart (const gsl_complex ex,
				   const gsl_complex ey,
				   const gsl_complex ez);

gsl_complex laser_polzn_vector_get (const laser_polzn_vector_t * e, const int p);

void laser_polzn_vector_set (laser_polzn_vector_t * e, const int p,
			     const gsl_complex val);

void laser_polzn_vector_dtor (laser_polzn_vector_t * e);

void laser_polzn_vector_rotate (laser_polzn_vector_t * e, const double phi,
				const double theta, const double chi);


/* Polarization tensor struct and functions. */

typedef struct _laser_polzn_tensor
{
  gsl_complex (*get) (const struct _laser_polzn_tensor * E, const int k,
		      const int p);
  void (*set) (struct _laser_polzn_tensor * E, const int k, const int p,
	       gsl_complex val);
  gsl_complex E[9];
} laser_polzn_tensor_t;

gsl_complex laser_polzn_tensor_get (const laser_polzn_tensor_t * E,
				    const int k, const int p);

void laser_polzn_tensor_set (laser_polzn_tensor_t * E, const int k, const int p,
			     const gsl_complex val);

laser_polzn_tensor_t *
laser_polzn_tensor_ctor_from_vectors (const laser_polzn_vector_t * e1,
				      const laser_polzn_vector_t * e2);

void laser_polzn_tensor_dtor (laser_polzn_tensor_t * E);

#endif /* __LASER_POLZN_H__ */
