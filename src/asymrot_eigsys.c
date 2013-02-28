#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "memory.h"
#include "jmarray_double.h"
#include "jkmarray_double.h"
#include "asymrot_eigsys.h"

struct _asymrot_eigsys
{
  int Jmax;
  double Bx, By, Bz;
  JKMarray_double_t *eigvec;
  JMarray_double_t *eigval;
};

static double
asymrot_eigsys_mtxel (const int J, const int K1, const int K2, 
			       const double Bx, const double By, 
			       const double Bz)
     /* Returns the matrix element of the asymmetric top Hamiltonian
        <J,K1|H|J,K2>.  Bx, By and Bz are the rotational constants associated
        with the three molecular axes - this is general, no assignation of the
        rotational constants A, B, C with the axes has been made here. See Zare
        page 272. */
{
  if (K1 == K2)
    {
      double a1 = J * (J + 1.0);
      double a2 = K1 * K1;
      return 0.5 * (a1 - a2) * (Bx + By) + Bz * a2;
    }

  else if (K2 == (K1 + 2))
    {
      double a1 = J * (J + 1.0);
      double a2 = K1 + 1.0;
      double a3 = K1 + 2.0;
      return 0.25 * (Bx - By) * sqrt ((a1 - K1 * a2) * (a1 - a2 * a3));
    }

  else if (K2 == (K1 - 2))
    {
      double a1 = J * (J + 1.0);
      double a2 = K1 - 1.0;
      double a3 = K1 - 2.0;
      return 0.25 * (Bx - By) * sqrt ((a1 - K1 * a2) * (a1 - a2 * a3));
    }

  else
    return 0.0;
}

static int
asymrot_eigsys_diagonalize (asymrot_eigsys_t * eigsys)
  /* Diagonalizes the asymmetric top Hamiltonian and stores the
     eigenvalues and eigenvector coefficients. */
{
  const double Bx = eigsys->Bx, By = eigsys->By, Bz = eigsys->Bz;
  const int Jmax = eigsys->Jmax;
  int J;

  for (J = 0; J <= Jmax; J++)
    {
      size_t dim = 2 * J + 1;
      gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc (dim);
      gsl_matrix *mtxin = gsl_matrix_calloc (dim, dim);	/* initially 0.0 */
      gsl_vector *eval = gsl_vector_alloc (dim);
      gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
      int status, Krow, Kcol, K, n;


      /* Calculate matrix elements (bra = row, ket = column). */
      for (Krow = -J; Krow <= J; Krow++)
	{
	  int irow = J + Krow; /* Adjust range to 0 .. 2J. */
	  
	  for (Kcol = -J; Kcol <= J; Kcol++)
	    {
	      int icol = J + Kcol; /* Adjust range to 0 .. 2J. */
	      double a = asymrot_eigsys_mtxel (J, Krow, Kcol, Bx, By, Bz);
	      gsl_matrix_set (mtxin, irow, icol, a);
	    }
	}

      /* Diagonalize */
      status = gsl_eigen_symmv (mtxin, eval, evec, workspace);
      if (status)
	{
	  fprintf (stderr,
		   "%s %d: error during diagonalization: %s\n",
		   __func__, __LINE__, gsl_strerror (status));
	  return -1;
	}

      /* Free workspace memory */
      gsl_eigen_symmv_free (workspace);
      gsl_matrix_free (mtxin);

      /* Sort according to eigenvalue - ascending order */
      status = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
      
      if (status)
	{
	  fprintf (stderr, "%s %d: error during sorting: %s\n",
		   __func__, __LINE__, gsl_strerror (status));
	  return -1;
	}
      
      /* Copy out eigenvalues and eigenvector coefficients. */
      for (n = -J; n <= J; n++)
	{
	  int nidx = J + n; /* Adjust range to 0 .. 2J. */
	  gsl_vector_view col = gsl_matrix_column (evec, nidx);
	  double val = gsl_vector_get (eval, nidx);

	  JMarray_double_set (eigsys->eigval, J, n, val);
	  
	  for (K = -J; K <= J; K++)
	    {
	      int Kidx = J + K; /* Adjust range to 0 .. 2J. */

	      val = gsl_vector_get (&col.vector, Kidx);
	      JKMarray_double_set (eigsys->eigvec, J, n, K, val);
	    }
	}
      
      /* Free local storage of eigenvalues and eigenvectors */
      gsl_vector_free (eval);
      gsl_matrix_free (evec);
    }

  return 0;
}

asymrot_eigsys_t *
asymrot_eigsys_ctor (const int Jmax, const double Bx, const double By, const double Bz)
{
  asymrot_eigsys_t * eigsys;

  if (MEMORY_ALLOC (eigsys) < 0)
    {
      MEMORY_OOMERR;
      fprintf (stderr, "%s %d: failed to allocate eigsys storage.\n",
	       __func__, __LINE__);
      return NULL;
    }
  
  eigsys->eigvec = JKMarray_double_ctor (Jmax);
  if (eigsys->eigvec == NULL)
    {
      MEMORY_OOMERR;
      fprintf (stderr, "%s %d: failed to allocate eigvec storage.\n",
	       __func__, __LINE__);
      MEMORY_FREE (eigsys);
      return NULL;
    }

  eigsys->eigval = JMarray_double_ctor (Jmax);
  if (eigsys->eigval == NULL)
    {
      MEMORY_OOMERR;
      fprintf (stderr, "%s %d: failed to allocate eigval storage.\n",
	       __func__, __LINE__);
      JKMarray_double_dtor (eigsys->eigvec);
      MEMORY_FREE (eigsys);
      return NULL;
    }

  eigsys->Jmax = Jmax;
  eigsys->Bx = Bx;
  eigsys->By = By;
  eigsys->Bz = Bz;

  if (asymrot_eigsys_diagonalize (eigsys) < 0)
    {
      MEMORY_OOMERR;
      fprintf (stderr, "%s %d: failed to diagonalize Hamiltonian.\n",
	       __func__, __LINE__);
      JKMarray_double_dtor (eigsys->eigvec);
      JMarray_double_dtor (eigsys->eigval);
      MEMORY_FREE (eigsys);
      return NULL;
    }

  return eigsys;
}

void
asymrot_eigsys_dtor (asymrot_eigsys_t *eigsys)
{
  JKMarray_double_dtor (eigsys->eigvec);
  JMarray_double_dtor (eigsys->eigval);
  MEMORY_FREE (eigsys);
}

double
asymrot_eigsys_eigvec_coef_get (asymrot_eigsys_t *eigsys, const int J, 
				const int n, const int K)
{
  return JKMarray_double_get (eigsys->eigvec, J, n, K);
}

double
asymrot_eigsys_eigval_get (asymrot_eigsys_t *eigsys, const int J, 
			   const int n)
{
  return JMarray_double_get (eigsys->eigval, J, n);
}

