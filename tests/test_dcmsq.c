#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "dcmsq.h"

int
main ()
{
  dcmsq_mtxel_t *mtxel;
  dcmsq_expval_t *expval1, *expval2;
  gsl_complex one = gsl_complex_rect (1.0, 0.0);

  mtxel = dcmsq_mtxel_ctor ();
  dcmsq_mtxel_calc (2, 2, 2, 2, 2, 2, mtxel);
  
  expval1 = dcmsq_expval_ctor ();
  expval2 = dcmsq_expval_ctor ();
  
  dcmsq_expval_add_mtxel_weighted_complex (expval1, mtxel, one);
  dcmsq_expval_zero (expval2);
  dcmsq_expval_add_weighted (expval2, expval1, 1.0);

  dcmsq_expval_dtor (expval1);
  dcmsq_expval_dtor (expval2);
  dcmsq_mtxel_dtor (mtxel);

  return 0;
}
