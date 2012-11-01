/* TODO: allow for 1/2 integer a.m. Use analytical formulas rather than 3js.*/

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include "dmtxel.h"


/* Useful warning messages. */
/* This function performs the normal angular momentum checks that would be
   performed in a sensible 3j calculation function as a safeguard. Because we're
   treating only integer angular momentum here, no test for 2(J+k+Jp) being even
   is required. */
/* #define DMTXEL_DEBUG 1 */

#ifdef DMTXEL_DEBUG
#include <stdio.h>
#define NAN 0.0/0.0 /* Not entirely portable. */
#endif /* DMTXEL_DEBUG */

double
dmtxel (const int J, const int K, const int M, const int Jp, const int Kp,
	const int Mp, const int k, const int p, const int q)
     /* Matrix element of the rotation matrix <J K M|D^k_pq|Jp Kp Mp>.  
	Currently only defined for integer angular momentum. */
{

#ifdef DMTXEL_DEBUG
  if (p != Mp - M)
    {
      fprintf (stderr, "dmtxel warning: p != Mp - M.\n");
      return 0.0;
    }
  else if (q != Kp - K)
    {
      fprintf (stderr, "dmtxel warning: q != Kp - K.\n");
      return 0.0;
    }
  else if (abs (q) > k)
    {
      fprintf (stderr, "dmtxel warning: abs(q) > k.\n");
      return NAN;
    }
  else if (abs (p) > k)
    {
      fprintf (stderr, "dmtxel warning: abs(p) > k.\n");
      return NAN;
    }
  else if (abs (M) > J)
    {
      fprintf (stderr, "dmtxel warning: abs(M) > J.\n");
      return NAN;
    }
  else if (abs (K) > J)
    {
      fprintf (stderr, "dmtxel warning: abs(K) > J.\n");
      return NAN;
    }
  else if (abs (Mp) > Jp)
    {
      fprintf (stderr, "dmtxel warning: abs(Mp) > Jp.\n");
      return NAN;
    }
  else if (abs (Kp) > Jp)
    {
      fprintf (stderr, "dmtxel warning: abs(Kp) > Jp.\n");
      return NAN;
    }
  else if ((Jp > J + k) || (Jp < abs (J - k)))
    {
      fprintf (stderr, "dmtxel warning: (J, k, Jp) not triangle.\n");
      return 0.0;
    }
#undef NAN
#endif /* DMTXEL_DEBUG */

  double a = sqrt ((2.0 * J + 1.0) * (2.0 * Jp + 1.0)) *
    gsl_sf_coupling_3j (2 * J, 2 * k, 2 * Jp, 2 * M, 2 * p, -2 * Mp) *
    gsl_sf_coupling_3j (2 * J, 2 * k, 2 * Jp, 2 * K, 2 * q, -2 * Kp);

  if (GSL_IS_ODD (Mp - Kp))	/* phase */
    return -a;
  else
    return a;
}
