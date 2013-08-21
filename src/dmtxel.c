/* TODO: Use analytical formulas rather than 3js.*/

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include "dmtxel.h"

#define __IS_ODD(n) ((n) & 1)

/* If DMTXEL_DEBUG is defined, this function performs the normal
   angular momentum checks that would be performed in a sensible 3j
   calculation function as a safeguard. Because we're treating only
   integer angular momentum here, no test for 2(J+k+Jp) being even is
   required. */
/* #define DMTXEL_DEBUG 1 */

#ifdef DMTXEL_DEBUG
#include <stdio.h>
#define __NAN 0.0/0.0		/* Not entirely portable. */
#endif /* DMTXEL_DEBUG */

double
dmtxel (const int two_J, const int two_K, const int two_M,
	const int two_Jp, const int two_Kp, const int two_Mp,
	const int two_k, const int two_p, const int two_q)
/* Matrix element of the rotation matrix <JKM|D^k_pq|JpKpMp>.
   Arguments passed are twice the angular momenta.*/
{

#ifdef DMTXEL_DEBUG
  if (two_p != two_Mp - two_M)
    {
      fprintf (stderr, "dmtxel warning: p != Mp - M.\n");
      return 0.0;
    }
  else if (two_q != two_Kp - two_K)
    {
      fprintf (stderr, "dmtxel warning: q != Kp - K.\n");
      return 0.0;
    }
  else if (abs (two_q) > two_k)
    {
      fprintf (stderr, "dmtxel warning: abs(q) > k.\n");
      return __NAN;
    }
  else if (abs (two_p) > two_k)
    {
      fprintf (stderr, "dmtxel warning: abs(p) > k.\n");
      return __NAN;
    }
  else if (abs (two_M) > two_J)
    {
      fprintf (stderr, "dmtxel warning: abs(M) > J.\n");
      return __NAN;
    }
  else if (abs (two_K) > two_J)
    {
      fprintf (stderr, "dmtxel warning: abs(K) > J.\n");
      return __NAN;
    }
  else if (abs (two_Mp) > two_Jp)
    {
      fprintf (stderr, "dmtxel warning: abs(Mp) > Jp.\n");
      return __NAN;
    }
  else if (abs (two_Kp) > two_Jp)
    {
      fprintf (stderr, "dmtxel warning: abs(Kp) > Jp.\n");
      return __NAN;
    }
  else if ((two_Jp > two_J + two_k) || (two_Jp < abs (two_J - two_k)))
    {
      fprintf (stderr, "dmtxel warning: (J, k, Jp) not triangle.\n");
      return 0.0;
    }
  else if (__IS_ODD (two_J + two_k + two_Jp))
    {
      fprintf (stderr, "dmtxel warning: 2(J+k+Jp) is odd.\n");
      return 0.0;
    }
#undef __NAN
#endif /* DMTXEL_DEBUG */

  double a = sqrt ((two_J + 1.0) * (two_Jp + 1.0)) *
    gsl_sf_coupling_3j (two_J, two_k, two_Jp, two_M, two_p, -two_Mp) *
    gsl_sf_coupling_3j (two_J, two_k, two_Jp, two_K, two_q, -two_Kp);

  if (__IS_ODD (two_Mp - two_Kp))	/* phase */
    return -a;
  else
    return a;
}
#undef __IS_ODD
