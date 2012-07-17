#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_coupling.h>

#include "laser_polzn.h"

void 
laser_polzn_vec_init (laser_polzn_vec *e)
{
  e->get = &laser_polzn_vec_get;
  e->set = &laser_polzn_vec_set;
  e->rotate = &laser_polzn_vec_rotate;
}


void 
laser_polzn_vec_free (laser_polzn_vec *e)
{
  e->get = NULL;
  e->set = NULL;
  e->rotate = NULL;
}


gsl_complex
laser_polzn_vec_get (const laser_polzn_vec *e, const int p)
{
  /* p takes the values (-1, 0, 1) stored in the array e[0, 1, 2]
     respectively */
  return e->e[p + 1];
}


void
laser_polzn_vec_set (laser_polzn_vec *e, const int p, const gsl_complex val)
{
  /* p takes the values (-1, 0, 1) stored in the array e[0, 1, 2]
     respectively */
  e->e[p + 1] = val;
}


void 
laser_polzn_vec_init_from_cart (laser_polzn_vec *e, const gsl_complex ex,
				const gsl_complex ey, const gsl_complex ez)
/* Takes the cartesian components of a complex vector and returns a complex
   vector containing the spherical tensor components. */
{
  const double rt2 = 1.0 / sqrt (2.0);
  gsl_complex i = gsl_complex_rect (0.0, 1.0);
  gsl_complex iey = gsl_complex_mul (ey, i);
  gsl_complex val;

  laser_polzn_vec_init (e);

  /* e_-1 =  1/sqrt(2) * (e_x - i*e_y) */
  val = gsl_complex_sub (ex, iey);
  val = gsl_complex_mul_real (val, rt2);
  e->set(e, -1, val);

  /* e_0 = e_z */
  e->set(e, 0, ez);

  /* e_+1 = - 1/sqrt(2) * (e_x + i*e_y) */
  val = gsl_complex_add (ex, iey);
  val = gsl_complex_mul_real (val, -rt2);
  e->set(e, 1, val);
}


static gsl_complex
laser_polzn_drot1 (const int p, const int q, const double phi, 
	     const double theta, const double chi)
{
  gsl_complex phiterm = gsl_complex_polar (1.0, -p * phi);
  gsl_complex chiterm = gsl_complex_polar (1.0, -q * chi);
  gsl_complex eterm = gsl_complex_mul (phiterm, chiterm);
  double d;
  // TODO: add checking of p and q values as being in range -1..1

  if (p == -1)
    {
      if (q == -1)
	d = 0.5 + 0.5 * cos(theta);
      else if (q == 0)
	d = sin(theta) / sqrt (2.0);
      else if (q == 1)
	d = 0.5 - 0.5 * cos(theta);
    }
  else if (p == 0) 
    {
      if (q == -1)
	d = -sin(theta) / sqrt (2.0);
      else if (q == 0)
	d = cos(theta);
      else if (q == 1)
	d = sin(theta) / sqrt (2.0);
    }
  else if (p == 1)
    {
      if (q == -1)
	d = 0.5 - 0.5 * cos(theta);
      else if (q == 0)
	d = -sin(theta) / sqrt (2.0);
      else if (q == 1)
	d = 0.5 + 0.5 * cos(theta);
    }
  
  return gsl_complex_mul_real (eterm, d);
      
}


void
laser_polzn_vec_rotate (laser_polzn_vec *e, const double phi,
			const double theta, const double chi)
/* Rotates the spherical components of the polarization vector by the euler
   angles phi, theta, chi. Spherical components of the unrotated vector of are
   passed in e, and the rotated components are returned in e. */
{
  laser_polzn_vec ein = *e;
  gsl_complex zero = gsl_complex_rect (0.0, 0.0);
  int p;
  
  /* Make a copy of the initial vector */
  ein = *e;

  for (p = -1; p <= 1; p++)
    {
      int q;

      e->set(e, p, zero);

      for (q = -1; q <= 1; q++)
	{
	  gsl_complex dmtx = laser_polzn_drot1 (p, q, phi, theta, chi);
	  gsl_complex term = gsl_complex_mul (dmtx, ein.get(&ein, q));
	  gsl_complex newval = gsl_complex_add (e->get(e, p), term); 
	  e->set(e, p, newval);
	}
    }
}


#define TIDX(k, p) ((k)*(k)+(k)+(p))

gsl_complex
laser_polzn_tensor_get (const laser_polzn_tensor *E, const int k, const int p)
{
  return E->E[TIDX (k, p)];
}

void
laser_polzn_tensor_set (laser_polzn_tensor *E, const int k, const int p,
			const gsl_complex val)
{
  E->E[TIDX (k, p)] = val;
}

#undef TIDX


void 
laser_polzn_tensor_init (laser_polzn_tensor *E)
{
  E->get = &laser_polzn_tensor_get;
  E->set = &laser_polzn_tensor_set;
}


void
laser_polzn_tensor_init_from_vecs (laser_polzn_tensor *E, const laser_polzn_vec *e1,
				   const laser_polzn_vec *e2)
/* Calculates the polarization tensor E = [e1* x e2]^k_q. Pass e1 to this
   function and not e1*. */
{
  int k;
  gsl_complex zero = gsl_complex_rect (0.0, 0.0);

  laser_polzn_tensor_init (E);

  for (k = 0; k <= 2; k++)
    {
      double a = sqrt (2.0 * k + 1.0);
      int q;

      for (q = -k; q <= k; q++)
	{
	  int p1;
	  double b;

	  if (GSL_IS_ODD (q))
	    b = -a;
	  else
	    b = a;

	  E->set(E, k, q, zero);

	  for (p1 = -1; p1 <= 1; p1++)
	    {
	      int p2 = q - p1;
	      double c;
	      gsl_complex d, e1cc;

	      if (abs (p2) > 1)
		continue;

	      c = b * gsl_sf_coupling_3j (2, 2, 2 * k, 2 * p1, 2 * p2, -2 * q);

	      e1cc = gsl_complex_conjugate (e1->get(e1, p1));

	      d = gsl_complex_mul (e1cc, e2->get(e2, p2));
	      d = gsl_complex_mul_real (d, c);

	      E->set(E, k, q, gsl_complex_add (E->get(E, k, q), d));
	    }
	}
    }
}

