#ifndef __LASER_POLZN_H__
#define __LASER_POLZN_H__ 1

#include <gsl/gsl_complex.h>

typedef struct _laser_polzn_vec
{
  gsl_complex (*get) (const struct _laser_polzn_vec * e, const int p);
  void (*set) (struct _laser_polzn_vec * e, const int p, gsl_complex val);
  void (*rotate) (struct _laser_polzn_vec * e,
		  const double phi, const double theta, const double chi);
  gsl_complex e[3];
} laser_polzn_vec_t;

void laser_polzn_vec_init_from_cart (laser_polzn_vec_t * e,
				     const gsl_complex ex,
				     const gsl_complex ey,
				     const gsl_complex ez);

gsl_complex laser_polzn_vec_get (const laser_polzn_vec_t * e, const int p);

void laser_polzn_vec_set (laser_polzn_vec_t * e, const int p,
			  const gsl_complex val);

void laser_polzn_vec_free (laser_polzn_vec_t * e);

void laser_polzn_vec_rotate (laser_polzn_vec_t * e, const double phi,
			     const double theta, const double chi);

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

void laser_polzn_tensor_init_from_vecs (laser_polzn_tensor_t * E,
					const laser_polzn_vec_t * e1,
					const laser_polzn_vec_t * e2);

void laser_polzn_tensor_free (laser_polzn_tensor_t * E);

#endif /* __LASER_POLZN_H__ */
