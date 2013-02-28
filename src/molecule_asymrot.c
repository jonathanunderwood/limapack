#include <stdlib.h>
#include <math.h>

#include "asymrot_eigsys.h"

typedef struct _asymrot_molecule
{
  double Bx, By, Bz;
  double partfn;
  double kT;
  int Jmax;
  asymrot_eigsys_t * eigsys;
  
} asymrot_molecule_t;

static void
molecule_asymrot_calc_partfn (asymrot_molecule_t * mol)
/* Calculate the rotational partition function. */
{
  int J;
  
  mol->partfn = 0.0;

  for (J = 0; J <= mol->Jmax; J++)
    {
      double dim = 2.0 * J + 1.0;
      int n;
      for (n = -J; n <= J; n++)
        {
          double energy = asymrot_eigsys_eigval_get (mol->eigsys, J, n);
          mol->partfn += dim * exp (-energy / mol->kT);
        }
    }

  return;
}
