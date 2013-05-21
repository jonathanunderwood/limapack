#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>

#include <libconfig.h>
#include <hdf5.h>

#ifdef BUILD_WITH_MPI
#include <mpi.h>
#endif

#include "molecule.h"
#include "laser.h"
#include "polarizability.h"
#include "jkmarray.h"
#include "jkmarray_int.h"
#include "memory.h"
#include "dmtxel.h"
#include "dcmsq.h"
#include "au.h"
#include "molecule_symrot.h"

struct _symrot_molecule
{
  molecule_t parent;
  double B_perp, B_par;
  double partfn;
  double kT;
  double coefmin;
  double poptol;
  int Jmax;
  int ncoef;
  int nexpval;
  polarizability_t *alpha;
  JKMarray_int_t *job_status;
};

typedef struct _symrot_molecule_tdse_worker
{
  molecule_tdse_worker_t parent;
  int J, M, K;
} symrot_molecule_tdse_worker_t;

typedef struct _symrot_molecule_expval
{
  dcmsq_expval_t *dcmsq;
} symrot_molecule_expval_t;

/* Note: in several of the functions below we don't make use of the
   first argument at present.  However, we keep the argument for
   generality - we may export these functions in the future to a
   library, or use them as call backs functions in other objects. */

static molecule_expval_t *
symrot_molecule_expval_ctor (const molecule_t * molecule)
{
  symrot_molecule_expval_t *expval = NULL;

  if (MEMORY_ALLOC (expval) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  expval->dcmsq = dcmsq_expval_ctor ();
  if (expval->dcmsq == NULL)
    {
      MEMORY_FREE (expval);
      return NULL;
    }

  return (molecule_expval_t *) expval;
}

static void
symrot_molecule_expval_dtor (const molecule_t * molecule,
			     molecule_expval_t * expval)
{
  symrot_molecule_expval_t *exp = (symrot_molecule_expval_t *) expval;

  dcmsq_expval_dtor (exp->dcmsq);
  MEMORY_FREE (exp);
}


static inline double
symrot_molecule_energy (const symrot_molecule_t * mol, const int J,
			const int K)
{
  return mol->B_perp * J * (J + 1) + (mol->B_par - mol->B_perp) * K * K;
}

static double
symrot_molecule_boltzmann_statwt (const symrot_molecule_t * mol, 
				  const int J, const int K)
/* Returns the statistical weight of a single (J, K, M) state in a
   Boltzmann distribution, hence no (2J+1) degeneracy pre-factor. */
{
  double E = symrot_molecule_energy (mol, J, K);

  return exp (-E / mol->kT) / mol->partfn;
}

static int
symrot_molecule_expval_fwrite (const molecule_t * molecule,
			       const molecule_expval_t * expval,
			       const hid_t * location)
{
  symrot_molecule_expval_t *exp = (symrot_molecule_expval_t *) expval;
  int ret = dcmsq_fwrite (exp->dcmsq, location);

  return ret;
}

static void
symrot_molecule_expval_zero (const symrot_molecule_t * molecule,
			     symrot_molecule_expval_t * expval)
/* Set all expectation values to zero. */
{
  dcmsq_expval_zero (expval->dcmsq);
}

static int
symrot_molecule_expval_calc (const molecule_t * molecule, const double *coef,
			     const double t, molecule_expval_t * expvalue)
/* Callback function used to calculate all needed expectation values. */
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;
  symrot_molecule_expval_t *expval = (symrot_molecule_expval_t *) expvalue;
  const double *coef_r = coef;
  const double *coef_i = coef + mol->ncoef;
  int J, Jmax = mol->Jmax;
  dcmsq_mtxel_t *dcmsq_mtxel = dcmsq_mtxel_ctor ();

  if (dcmsq_mtxel == NULL)
    return -1;

  /* Initialize all expectation values to zero. */
  symrot_molecule_expval_zero (mol, expval);

  for (J = 0; J <= Jmax; J++)
    {
      int Jpmin = J - 2;
      int Jpmax;
      int M;

      if (Jpmin < 0)
	Jpmin = 0;

      Jpmax = J + 2;
      if (Jpmax > Jmax)
	Jpmax = Jmax;

      for (M = -J; M <= J; M++)
	{
	  int K;

	  for (K = -J; K <= J; K++)
	    {

	      double E = symrot_molecule_energy (mol, J, K);
	      int i1 = JKMarray_idx (J, K, M);
	      int Jp;
	      gsl_complex A;

	      GSL_SET_COMPLEX (&A, coef_r[i1], coef_i[i1]);

	      if (gsl_complex_abs (A) < mol->coefmin)
		continue;

	      A = gsl_complex_conjugate (A);

	      for (Jp = Jpmin; Jp <= Jpmax; Jp++)
		{
		  int p;

		  for (p = -2; p <= 2; p++)
		    {
		      int Mp = M + p;
		      int q;

		      if (abs (Mp) > Jp)
			continue;

		      for (q = -2; q <= 2; q += 2)
			{
			  int Kp = K + q;
			  double Ep;
			  int i2;
			  gsl_complex Ap, c1, c2, c3;

			  if (abs (Kp) > Jp)
			    continue;

			  i2 = JKMarray_idx (Jp, Kp, Mp);

			  GSL_SET_COMPLEX (&Ap, coef_r[i2], coef_i[i2]);

			  if (gsl_complex_abs (Ap) < mol->coefmin)
			    continue;

			  Ep = symrot_molecule_energy (mol, Jp, Kp);

			  c1 = gsl_complex_mul (A, Ap);
			  c2 = gsl_complex_polar (1.0, (E - Ep) * t);
			  c3 = gsl_complex_mul (c1, c2);

			  dcmsq_mtxel_calc (dcmsq_mtxel, J, K, M, Jp, Kp, Mp);

			  dcmsq_expval_add_mtxel_weighted_complex
			    (expval->dcmsq, dcmsq_mtxel, c3);
			}
		    }
		}
	    }
	}
    }

  dcmsq_mtxel_dtor (dcmsq_mtxel);

  return 0;
}

static void
symrot_molecule_expval_add_weighted (const molecule_t * molecule,
				     molecule_expval_t * a,
				     const molecule_expval_t * b,
				     const double weight)
/* Perform a = a + (weight * b) */
{
  symrot_molecule_expval_t *aa = (symrot_molecule_expval_t *) a;
  symrot_molecule_expval_t *bb = (symrot_molecule_expval_t *) b;

  dcmsq_expval_add_weighted (aa->dcmsq, bb->dcmsq, weight);
  return;
}

#ifdef BUILD_WITH_MPI
static int
symrot_molecule_expval_mpi_send (const molecule_t * molecule,
				 const molecule_expval_t * expval,
				 int dest, int tag, MPI_Comm comm)
{
  const symrot_molecule_expval_t *ev = (symrot_molecule_expval_t *) expval;
  int ret;

  ret = dcmsq_expval_mpi_send (ev->dcmsq, dest, tag, comm);

  return ret;
}

static int
symrot_molecule_expval_mpi_recv (const molecule_t * molecule,
				 molecule_expval_t * expval,
				 int dest, int tag, MPI_Comm comm)
{
  symrot_molecule_expval_t *ev = (symrot_molecule_expval_t *) expval;
  int ret;

  ret = dcmsq_expval_mpi_recv (ev->dcmsq, dest, tag, comm);

  return ret;
}
#endif

/* These macros define TDSE job states for the following functions. */
#define __TODO 0
#define __STARTED 1
#define __DONE 2

static int
symrot_molecule_get_tdse_job (molecule_t * molecule,
			      molecule_tdse_worker_t * worker)
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;
  symrot_molecule_tdse_worker_t *w = (symrot_molecule_tdse_worker_t *) worker;
  int J, M, K;

  for (J = 0; J <= mol->Jmax; J++)
    for (K = -J; K <= J; K++)
      for (M = -J; M <= J; M++)
	{
	  int status = JKMarray_int_get (mol->job_status, J, K, M);
	  if (status == __TODO)
	    {
	      double wt = symrot_molecule_boltzmann_statwt (mol, J, K);

	      if (wt < mol->poptol)	/* Population in this state negligible */
		JKMarray_int_set (mol->job_status, J, K, M, __DONE);
	      else
		{
		  JKMarray_int_set (mol->job_status, J, K, M, __STARTED);
		  w->J = J;
		  w->K = K;
		  w->M = M;
		  snprintf (w->parent.description,
			    MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH,
			    "J: %d K: %d M: %d", J, K, M);
		  return 0;
		}
	    }
	}

  /* No jobs in the TODO state. */
  return -1;
}

void
symrot_molecule_set_tdse_job_done (molecule_t * molecule,
				   molecule_tdse_worker_t * worker)
{
  symrot_molecule_t *m = (symrot_molecule_t *) molecule;
  symrot_molecule_tdse_worker_t *w = (symrot_molecule_tdse_worker_t *) worker;

  JKMarray_int_set (m->job_status, w->J, w->K, w->M, __DONE);
}

static void
symrot_molecule_get_tdse_job_coef (const molecule_t * molecule,
				   const molecule_tdse_worker_t * worker,
				   double *coef)
{
  symrot_molecule_tdse_worker_t *w = (symrot_molecule_tdse_worker_t *) worker;
  int i;

  /* Set all real and imaginary parts of the wavefunction coefficients
     to zero. */
  for (i = 0; i < 2 * molecule->get_ncoef (molecule); i++)
    coef[i] = 0.0;

  /* And now set the real part of the (J, K, M) coefficient to 1. */
  coef[JKMarray_idx (w->J, w->K, w->M)] = 1.0;	/* Im part is 0.0, set above. */

  return;
}

static double
symrot_molecule_get_tdse_job_weight (const molecule_t * molecule,
				     const molecule_tdse_worker_t * worker)
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;
  symrot_molecule_tdse_worker_t *w = (symrot_molecule_tdse_worker_t *) worker;

  return symrot_molecule_boltzmann_statwt (mol, w->J, w->K);
}

molecule_tdse_worker_t *
symrot_molecule_tdse_worker_ctor (const molecule_t * molecule)
{
  symrot_molecule_tdse_worker_t *worker;

  if (MEMORY_ALLOC (worker) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  return (molecule_tdse_worker_t *) worker;
}

void
symrot_molecule_tdse_worker_dtor (const molecule_t * molecule,
				  molecule_tdse_worker_t * worker)
{
  MEMORY_FREE (worker);
}

#ifdef BUILD_WITH_MPI
static int
symrot_molecule_tdse_worker_mpi_send (const molecule_t * self,
				      const molecule_tdse_worker_t * worker,
				      int dest, int tag, MPI_Comm comm)
{
  const symrot_molecule_tdse_worker_t *w =
    (symrot_molecule_tdse_worker_t *) worker;
  int ret;

  /* The casts to (void *) here are because the MPI API doesn't
     include the const qualifier, because academic programmers are
     morons. So, we have to discard the const qualifier. This is
     apparently fixed in the MPI 3 spec. */
  ret = MPI_Send ((void *) &(w->parent.description),
		  MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH,
		  MPI_CHAR, dest, tag, comm);
  ret += MPI_Send ((void *) &(w->J), 1, MPI_INT, dest, tag, comm);
  ret += MPI_Send ((void *) &(w->K), 1, MPI_INT, dest, tag, comm);
  ret += MPI_Send ((void *) &(w->M), 1, MPI_INT, dest, tag, comm);

  return ret;
}

static int
symrot_molecule_tdse_worker_mpi_recv (const molecule_t * self,
				      molecule_tdse_worker_t * worker,
				      int dest, int tag, MPI_Comm comm)
{
  symrot_molecule_tdse_worker_t *w = (symrot_molecule_tdse_worker_t *) worker;
  int ret;

  ret = MPI_Recv ((void *) &(w->parent.description),
		  MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH,
		  MPI_CHAR, dest, tag, comm, MPI_STATUS_IGNORE);
  ret += MPI_Recv (&(w->J), 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);
  ret += MPI_Recv (&(w->K), 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);
  ret += MPI_Recv (&(w->M), 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);

  return ret;
}
#endif

int
symrot_molecule_get_ncoef (const molecule_t * molecule)
{
  symrot_molecule_t *m = (symrot_molecule_t *) molecule;

  return m->ncoef;
}

int
symrot_molecule_tdse_rhs (const molecule_t * molecule,
			  const laser_collection_t * lasers,
			  const double t, const double *coef, double *deriv)
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;
  const double *coef_real = coef, *coef_imag = coef + mol->ncoef;
  double *deriv_real = deriv, *deriv_imag = deriv + mol->ncoef;
  int Jmax = mol->Jmax;
  int ilas;
  laser_polzn_tensor_t *E = laser_polzn_tensor_ctor ();

  if (E == NULL)
    return GSL_FAILURE;

  for (ilas = 0; ilas < lasers->nlasers; ilas++)
    {
      int jlas;
      double env1 =
	lasers->laser[ilas]->get_envelope (lasers->laser[ilas], t);

      if (env1 < 0.0)
	continue;

      for (jlas = 0; jlas < lasers->nlasers; jlas++)
	{
	  double f, env2 =
	    lasers->laser[jlas]->get_envelope (lasers->laser[jlas], t);
	  int J;

	  if (env2 < 0.0)
	    continue;

	  f = env1 * env2;

	  laser_polzn_tensor_set_from_vectors
	    (E,
	     lasers->laser[ilas]->get_polzn_vector (lasers->laser[ilas], t),
	     lasers->laser[jlas]->get_polzn_vector (lasers->laser[jlas], t));

	  for (J = 0; J <= Jmax; J++)
	    {
	      int Jpmin = J - 2;
	      int Jpmax, M;

	      if (Jpmin < 0)
		Jpmin = 0;

	      Jpmax = J + 2;
	      if (Jpmax > Jmax)
		Jpmax = Jmax;

	      for (M = -J; M <= J; M++)
		{
		  int K;

		  for (K = -J; K <= J; K++)
		    {
		      int Jp;
		      int i1 = JKMarray_idx (J, K, M);
		      double E_JK = symrot_molecule_energy (mol, J, K);

		      for (Jp = Jpmin; Jp <= Jpmax; Jp++)
			{
			  int Mp, Mpmin = GSL_MAX (-Jp, M - 2);
			  int Mpmax = GSL_MIN (Jp, M + 2);
			  int Kpmin = -Jp;
			  int Kpmax = Jp;

			  for (Mp = Mpmin; Mp <= Mpmax; Mp++)
			    {
			      int Kp;

			      for (Kp = Kpmin; Kp <= Kpmax; Kp++)
				{
				  int i2 = JKMarray_idx (Jp, Kp, Mp);
				  double E_JpKp;
				  gsl_complex Ap, c1, eterm;
				  int k, kmin, p = Mp - M;

				  GSL_SET_COMPLEX
				    (&Ap, coef_real[i2], coef_imag[i2]);

				  if (gsl_complex_abs (Ap) < mol->coefmin)
				    continue;

				  E_JpKp =
				    symrot_molecule_energy (mol, Jp, Kp);

				  eterm =
				    gsl_complex_polar (1.0,
						       (E_JK - E_JpKp) * t);

				  c1 = gsl_complex_mul (Ap, eterm);

				  if ((J == Jp) && (p == 0))
				    kmin = 0;
				  else
				    kmin = 2;

				  for (k = kmin; k <= 2; k += 2)
				    {
				      gsl_complex Ekp = E->get (E, k, p);
				      gsl_complex c2 =
					gsl_complex_mul (c1, Ekp);
				      int q;

				      for (q = -k; q <= k; q += 2)
					{
					  int Kp = K + q;
					  double b, alpha;
					  gsl_complex c3;

					  if (abs (Kp) > Jp)
					    continue;

					  alpha =
					    polarizability_get (mol->alpha, k,
								-q);

					  b = -alpha * f *
					    dmtxel (J, K, M, Jp, Kp, Mp,
						    k, p, q);

					  c3 = gsl_complex_mul_real (c2, b);

					  /* Signs here reflect the factor of i in the 
					     TDSE. */
					  deriv_real[i1] += GSL_IMAG (c3);
					  deriv_imag[i1] -= GSL_REAL (c3);
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

  laser_polzn_tensor_dtor (E);

  return GSL_SUCCESS;
}

int
symrot_molecule_check_populations (const molecule_t * molecule,
				   const double *coef)
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;
  const double *coef_r = coef;
  const double *coef_i = coef + mol->ncoef;
  int J;

  /* Check populations in Jmax manifold, and also Jmax - 1
     manifold. This is actually crude, as it's probable that the
     highest energy state in the field at any given point may not
     correspond tto those of Jmax - n dependence is unclear. Something
     to fix for the future. */
  for (J = mol->Jmax - 1; J <= mol->Jmax; J++)
    {
      int K;

      for (K = -J; K <= J; K++)
	{
	  int M;
	  double re, im, pop = 0.0;

	  for (M = -J; M <= J; M++)
	    {
	      int idx = JKMarray_idx (J, K, M);
	      re = coef_r[idx];
	      im = coef_i[idx];
	      pop += (re * re) + (im * im);	/* Population in this (J, M) state */
	    }
	  /* If total population in this J manifold is greater than
	     poptol, then it's considered unacceptable. */
	  if (pop > mol->poptol)
	    {
	      fprintf (stderr,
		       "Population greater than poptol in J=%d. Population=%g, poptol=%g\n",
		       J, pop, mol->poptol);
	      return -1;
	    }
	}
    }

  return 0;
}

symrot_molecule_t *
symrot_molecule_ctor (const double B_par, const double B_perp,
		      const int Jmax, const double T,
		      const double alpha_par, const double alpha_perp,
		      const double coefmin, const double poptol)
/* B_par, B_perp are rotational constants in wavenumbers associated
   with rotation about axes parallel and perpendicular to the molecule
   z axis respectively.
   
   Jmax is the maximum value of J to consider.

   alpha_par, alpha_perp are the molecular frame components of
   polarizability parallalel and perpendicular to the molecule z-axis
   in Angstrom cubed.
   
   coefmin is the minimum value for a eigenstate coefficient which we
   bother to propagate in the TDSE.

   poptol is the maximum allowed population in the Jmax and Jmax-1
   manifolds during propagation. If the population in these manifolds
   exceeds poptol the calculation will bail out.
*/
{
  double a_par, a_perp;
  symrot_molecule_t *mol;
  int J;

  if (MEMORY_ALLOC (mol) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  /* Allocate parent structure and register dispatch functions. */
  molecule_dispatch_register ((molecule_t *) mol,
			      symrot_molecule_tdse_rhs,
			      symrot_molecule_get_tdse_job,
			      symrot_molecule_get_tdse_job_coef,
			      symrot_molecule_get_tdse_job_weight,
			      symrot_molecule_set_tdse_job_done,
			      symrot_molecule_check_populations,
			      symrot_molecule_tdse_worker_ctor,
			      symrot_molecule_tdse_worker_dtor,
			      symrot_molecule_get_ncoef,
			      symrot_molecule_expval_ctor,
			      symrot_molecule_expval_dtor,
			      symrot_molecule_expval_calc,
			      symrot_molecule_expval_add_weighted,
			      symrot_molecule_expval_fwrite,
#ifdef BUILD_WITH_MPI
			      symrot_molecule_expval_mpi_send,
			      symrot_molecule_expval_mpi_recv,
			      symrot_molecule_tdse_worker_mpi_send,
			      symrot_molecule_tdse_worker_mpi_recv,
#endif
			      symrot_molecule_dtor);
  mol->Jmax = Jmax;
  mol->poptol = poptol;
  mol->coefmin = coefmin;

  /* Convert to atomic units. */
  mol->B_par = WN_TO_AU (B_par);
  mol->B_perp = WN_TO_AU (B_perp);


  /* Caulculate kT in atomic units from T. */
  mol->kT = T_TO_KTAU (T);

  /* Calculate polarizability tensor in atomic units. */
  a_par = ANG3_TO_AU (alpha_par);
  a_perp = ANG3_TO_AU (alpha_perp);
  mol->alpha = polarizability_from_cart_ctor (a_perp, a_perp, a_par);

  if (mol->alpha == NULL)
    {
      fprintf (stderr, "%s %d: failed to allocate memory for mol->alpha.\n",
	       __func__, __LINE__);
      symrot_molecule_dtor ((molecule_t *) mol);
      return NULL;
    }

  /* Number of complex coefficients. */
  mol->ncoef = JKMarray_dim (mol->Jmax);

  /* Calculate and store the partition function. */
  mol->partfn = 0.0;
  for (J = 0; J <= mol->Jmax; J++)
    {
      int K;

      for (K = -J; K <= J; K++)
	{
	  double dim = 2.0 * J + 1.0;
	  double E = symrot_molecule_energy (mol, J, K);
	  mol->partfn += dim * exp (-E / mol->kT);
	}
    }

  /* Setup array for storing TDSE job status. */
  mol->job_status = JKMarray_int_ctor (mol->Jmax);
  if (mol->job_status == NULL)
    {
      fprintf (stderr,
	       "%s %d: failed to allocate memory for mol->tdse_status.\n",
	       __func__, __LINE__);
      symrot_molecule_dtor ((molecule_t *) mol);
      return NULL;
    }

  for (J = 0; J <= mol->Jmax; J++)
    {
      int K;

      for (K = -J; K <= J; K++)
	{
	  int M;

	  for (M = -J; M <= J; M++)
	    JKMarray_int_set (mol->job_status, J, K, M, __TODO);
	}
    }

  return mol;
}

molecule_t *
symrot_molecule_cfg_parse_ctor (const config_setting_t * cfg)
{
  symrot_molecule_t *mol;
  double B_par, B_perp;
  double alpha_par, alpha_perp;
  double T, coefmin, poptol;
  int Jmax;

  if (!(config_setting_lookup_float (cfg, "B_par", &B_par) &&
	config_setting_lookup_float (cfg, "B_perp", &B_perp) &&
	config_setting_lookup_float (cfg, "alpha_par", &alpha_par) &&
	config_setting_lookup_float (cfg, "alpha_perp", &alpha_perp) &&
	config_setting_lookup_float (cfg, "T", &T) &&
	config_setting_lookup_float (cfg, "coefmin", &coefmin) &&
	config_setting_lookup_float (cfg, "poptol", &poptol) &&
	config_setting_lookup_int (cfg, "Jmax", &Jmax)))
    {
      fprintf (stderr, "%s %d: incomplete symrot molecule configuration.\n",
	       __func__, __LINE__);
      return NULL;
    }

  mol = symrot_molecule_ctor (B_par, B_perp, Jmax, T, alpha_par, alpha_perp,
			      coefmin, poptol);

  if (mol == NULL)
    fprintf (stderr,
	     "%s %d: failed to create symrot_molecule structure after parsing config.\n",
	     __func__, __LINE__);

  return (molecule_t *) mol;
}

void
symrot_molecule_dtor (molecule_t * molecule)
{
  symrot_molecule_t *mol = (symrot_molecule_t *) molecule;

  if (mol->job_status != NULL)
    JKMarray_int_dtor (mol->job_status);

  if (mol->alpha != NULL)
    polarizability_dtor (mol->alpha);

  MEMORY_FREE (mol);
}

#undef __TODO
#undef __STARTED
#undef __DONE
