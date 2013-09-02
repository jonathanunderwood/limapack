#include <stdlib.h>
#include <stdio.h>
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
#include "jmarray.h"
#include "jmarray_int.h"
#include "memory.h"
#include "molecule_linear.h"
#include "dmtxel.h"
#include "dcmsq.h"
#include "au.h"

struct _linear_molecule
{
  molecule_t parent;
  double B;
  double kT, partfn;
  double coefmin;
  double poptol;
  double oddwt, evenwt;
  polarizability_t * alpha;
  int two_Jmax;
  int ncoef;
  JMarray_int_t * job_status;
};

typedef struct _linear_molecule_tdse_worker
{
  molecule_tdse_worker_t parent;
  int two_J, two_M;
} linear_molecule_tdse_worker_t;

typedef struct _linear_molecule_expval
{
  dcmsq_expval_t *dcmsq;
} linear_molecule_expval_t;


static molecule_expval_t *
linear_molecule_expval_ctor (const molecule_t *molecule)
/* Note at present first argument is unused. Kept for generality and future use. */
{
  linear_molecule_expval_t *expval = NULL;

  if (MEMORY_ALLOC(expval) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  expval->dcmsq = dcmsq_expval_ctor();
  if (expval->dcmsq == NULL)
    {
      MEMORY_FREE(expval);
      return NULL;
    }

  return (molecule_expval_t *) expval;
}

static void 
linear_molecule_expval_dtor (const molecule_t *molecule, molecule_expval_t *expval)
/* Note at present first argument is unused. Kept for generality and future use. */
{
  linear_molecule_expval_t *exp = (linear_molecule_expval_t *) expval;

  dcmsq_expval_dtor (exp->dcmsq);
  MEMORY_FREE(exp);
}


static inline double
linear_molecule_energy(const linear_molecule_t *mol, const int two_J)
{
  double J = 0.5 * two_J;
  return mol->B * J * (J + 1);
}

static double
linear_molecule_boltzmann_statwt(const linear_molecule_t * mol,
				 const int two_J)
/* Returns the statistical weight of a single (J,M) state in a
   Boltzmann distribution, hence no (2J+1) degeneracy pre-factor. 
   FIXME: At present this deals only with integer angular momenta.
*/
{ 
  double wt = GSL_IS_ODD (two_J / 2) ? mol->oddwt : mol->evenwt; 
  double E = linear_molecule_energy (mol, two_J); 

  return wt * exp (-E / mol->kT) / mol->partfn; 
}

static int
linear_molecule_expval_fwrite (const molecule_t *molecule,
			       const molecule_expval_t *expval,
			       const hid_t *location)
{
  // TODO: want a metadata node that contains the config?
  linear_molecule_expval_t *exp = (linear_molecule_expval_t *) expval;
  int ret = dcmsq_fwrite (exp->dcmsq, location);
  return ret;
}

static void 
linear_molecule_expval_zero(const linear_molecule_t *molecule, 
			    linear_molecule_expval_t *expval)
/* Set all expectation values to zero. */
{
  dcmsq_expval_zero(expval->dcmsq);
}

static int
linear_molecule_expval_calc (const molecule_t *molecule, const double *coef, 
			     const double t, molecule_expval_t *expvalue)
/* Callback function used to calculate all needed expectation values. */
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  linear_molecule_expval_t *expval = (linear_molecule_expval_t *) expvalue;
  const double *coef_r = coef;
  const double *coef_i = coef + mol->ncoef;
  int J;
  dcmsq_mtxel_t *dcmsq_mtxel = dcmsq_mtxel_ctor();

  if (dcmsq_mtxel == NULL)
    return -1;

  /* Initialize all expectation values to zero. */
  linear_molecule_expval_zero(mol, expval);

  for (two_J = 0; two_J <= mol->two_Jmax; two_J += 2)
    {
      int two_Jpmin = two_J - 4, two_Jpmax = two_J + 4, two_M;

      if (two_Jpmin < 0)
	two_Jpmin = 0;

      if (two_Jpmax > mol->two_Jmax)
	two_Jpmax = mol->two_Jmax;

      for (two_M = -two_J; two_M <= two_J; two_M += 2)
	{
	  double E = linear_molecule_energy (mol, two_J);
	  int i1 = JMarray_idx (two_J, two_M);
	  int two_Jp;
	  gsl_complex A;

	  GSL_SET_COMPLEX (&A, coef_r[i1], coef_i[i1]);

	  if (gsl_complex_abs (A) < mol->coefmin)
	    continue;

	  A = gsl_complex_conjugate (A);

	  for (two_Jp = two_Jpmin; two_Jp <= two_Jpmax; Jp += 2)
	    {
	      int two_p;

	      for (two_p = -4; two_p <= 4; two_p += 2)
		{
		  int two_Mp = two_M + two_p;
		  double Ep;
		  int i2;
		  gsl_complex Ap, c1, c2, c3;


		  if (abs (two_Mp) > two_Jp)
		    continue;

		  Ep = linear_molecule_energy (mol, two_Jp);
		  i2 = JMarray_idx (two_Jp, two_Mp);

		  GSL_SET_COMPLEX (&Ap, coef_r[i2], coef_i[i2]);

		  if (gsl_complex_abs (Ap) < mol->coefmin)
		    continue;

		  c1 = gsl_complex_mul (A, Ap);
		  c2 = gsl_complex_polar (1.0, (E - Ep) * t);
		  c3 = gsl_complex_mul (c1, c2);

		  dcmsq_mtxel_calc (dcmsq_mtxel, two_J, 0, two_M, two_Jp, 0, two_Mp);
		  dcmsq_expval_add_mtxel_weighted_complex (expval->dcmsq, dcmsq_mtxel, c3);
		}
	    }
	}
    }

  dcmsq_mtxel_dtor(dcmsq_mtxel);

  return 0;
}

static void 
linear_molecule_expval_add_weighted (const molecule_t *molecule, 
				     molecule_expval_t *a, const molecule_expval_t *b,
				     const double weight)
/* Perform a = a + (weight * b). At present first argument is unused,
   but kept for future use. */ 
{
  linear_molecule_expval_t *aa = (linear_molecule_expval_t *) a;
  linear_molecule_expval_t *bb = (linear_molecule_expval_t *) b;
  
  dcmsq_expval_add_weighted (aa->dcmsq, bb->dcmsq, weight);
  return;
}

#ifdef BUILD_WITH_MPI
static int
linear_molecule_expval_mpi_send(const molecule_t *molecule,
				const molecule_expval_t *expval, 
				int dest, int tag, MPI_Comm comm)
/* First argument presently unused. */
{
  const linear_molecule_expval_t *ev = (linear_molecule_expval_t *) expval;
  int ret;

  ret = dcmsq_expval_mpi_send(ev->dcmsq, dest, tag, comm);

  return ret;
}

static int
linear_molecule_expval_mpi_recv(const molecule_t *molecule, 
				molecule_expval_t *expval, 
				int dest, int tag, MPI_Comm comm)
  /* First argument presently unused. */
{
  linear_molecule_expval_t *ev = (linear_molecule_expval_t *) expval;
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
linear_molecule_get_tdse_job(molecule_t *molecule, 
			     molecule_tdse_worker_t *worker)
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;
  int two_J, two_M;

  for (two_J = 0; two_J <= mol->two_Jmax; J += 2)
    for (two_M = -two_J; two_M <= two_J; two_M += 2)
      {
	int status = 
	  JMarray_int_get(mol->job_status, two_J, two_M);
	if (status == __TODO)
	  {
	    double wt = linear_molecule_boltzmann_statwt(mol, two_J);

	    if (wt < mol->poptol) /* Population in this state negligible */
		JMarray_int_set(mol->job_status, two_J, two_M, __DONE);
	    else
	      {
		JMarray_int_set(mol->job_status, two_J, two_M, __STARTED);
		w->two_J = two_J;
		w->two_M = two_M;
		snprintf(w->parent.description, 
			 MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH, 
			 "J: %.1f M: %.1f", two_J * 0.5, two_M * 0.5);
		return 0;
	      }
	  }
      }

  /* No jobs in the TODO state. */
  return -1;
}

void
linear_molecule_set_tdse_job_done (molecule_t *molecule, 
				   molecule_tdse_worker_t *worker)
{
  linear_molecule_t *m = 
    (linear_molecule_t *) molecule;
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;

  JMarray_int_set(m->job_status, w->two_J, w->two_M, __DONE);
}

static void
linear_molecule_get_tdse_job_coef(const molecule_t *molecule, 
				  const molecule_tdse_worker_t *worker,
				  double *coef)
{
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;
  int i;

  /* Set all real and imaginary parts of the wavefunction coefficients
     to zero. */ 
  for (i = 0; i < 2 * molecule->get_ncoef(molecule); i++)
    coef[i] = 0.0;
		
  /* And now set the real part of the (J, M) coefficient to 1. */
  coef[JMarray_idx (w->two_J, w->two_M)] = 1.0; /* Im part is 0.0, set above. */

  return;
}

static double
linear_molecule_get_tdse_job_weight(const molecule_t *molecule, 
				    const molecule_tdse_worker_t *worker)
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;
  
  return linear_molecule_boltzmann_statwt(mol, w->two_J);
}



molecule_tdse_worker_t * 
linear_molecule_tdse_worker_ctor(const molecule_t *molecule)
/* Note: at present the passed argument is not used. */
{
  linear_molecule_tdse_worker_t * worker;

  if (MEMORY_ALLOC(worker) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  return (molecule_tdse_worker_t *) worker;
}

void 
linear_molecule_tdse_worker_dtor(const molecule_t *molecule, 
				 molecule_tdse_worker_t *worker)
/* Note at present the first argument is not used. */
{
  /* molecule_tdse_worker_dtor(worker->parent); */
  MEMORY_FREE(worker);
}


#ifdef BUILD_WITH_MPI
static int 
linear_molecule_tdse_worker_mpi_send (const molecule_t *self, 
				      const molecule_tdse_worker_t *worker,
				      int dest, int tag, MPI_Comm comm)
{
  const linear_molecule_tdse_worker_t *w =
    (linear_molecule_tdse_worker_t *) worker;
  int ret;

  /* The casts to (void *) here are because the MPI API doesn't
     include the const qualifier, because academic programmers are
     morons. So, we have to discard the const qualifier. This is
     apparently fixed in the MPI 3 spec. */
  ret = MPI_Send ((void *) &(w->parent.description), 
		  MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH, 
		  MPI_CHAR, dest, tag, comm);
  ret += MPI_Send ((void *) &(w->two_J), 1, MPI_INT, dest, tag, comm);
  ret += MPI_Send ((void *) &(w->two_M), 1, MPI_INT, dest, tag, comm);

  return ret;
}

static int 
linear_molecule_tdse_worker_mpi_recv (const molecule_t *self, 
				      molecule_tdse_worker_t *worker,
				      int dest, int tag, MPI_Comm comm)
{
  linear_molecule_tdse_worker_t *w =
    (linear_molecule_tdse_worker_t *) worker;
  int ret;

  ret = MPI_Recv ((void *) &(w->parent.description), 
		  MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH, 
		  MPI_CHAR, dest, tag, comm, MPI_STATUS_IGNORE);
  ret += MPI_Recv (&(w->two_J), 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);
  ret += MPI_Recv (&(w->two_M), 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);

  return ret;
}
#endif

int
linear_molecule_get_ncoef(const molecule_t * molecule)
{
  linear_molecule_t * m = (linear_molecule_t *) molecule;

  return m->ncoef;
}

int 
linear_molecule_tdse_rhs(const molecule_t *molecule,
			 const laser_collection_t* lasers,
			 const double t,
			 const double *coef,
			 double *deriv)
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  const double *coef_real = coef, *coef_imag = coef + mol->ncoef;
  double *deriv_real = deriv, *deriv_imag = deriv + mol->ncoef;
  int ilas;
  laser_polzn_tensor_t * E = laser_polzn_tensor_ctor();

  if (E == NULL)
    return GSL_FAILURE;

  for (ilas = 0; ilas < lasers->nlasers; ilas++)
    {
      int jlas;
      double env1 = lasers->laser[ilas]->get_envelope(lasers->laser[ilas], t);
      
      if (env1 < 0.0)
	continue;

      for (jlas = 0; jlas < lasers->nlasers; jlas++)
	{
	  double f, env2 = lasers->laser[jlas]->get_envelope(lasers->laser[jlas], t);
	  int two_J;//, idx1 = -1;

	  if (env2 < 0.0)
	    continue;

	  f = env1 * env2;

	  laser_polzn_tensor_set_from_vectors
	    (E, lasers->laser[ilas]->get_polzn_vector(lasers->laser[ilas], t),
	     lasers->laser[jlas]->get_polzn_vector(lasers->laser[jlas], t));

	  for (two_J = 0; two_J <= mol->two_Jmax; two_J += 2)
	    {
	      int two_Jpmin = two_J - 4, two_Jpmax = two_J + 4, two_M;
	      double EJ;

	      if (two_Jpmin < 0)
		two_Jpmin = 0;

	      if (two_Jpmax > mol->two_Jmax)
		two_Jpmax = mol->two_Jmax;

	      EJ = linear_molecule_energy (mol, two_J);
	      
	      for (two_M = -two_J; two_M <= two_J; two_M += 2)
		{
		  int Jp;
		  int i1 = JMarray_idx (two_J, two_M);
		  
		  for (two_Jp = two_Jpmin; two_Jp <= two_Jpmax; two_Jp += 2)
		    {
		      int Mp, Mpmin = GSL_MAX (-two_Jp, two_M - 4);
		      int Mpmax = GSL_MIN (two_Jp, two_M + 4);
		      double EJp = linear_molecule_energy (mol, two_Jp);

		      for (two_Mp = two_Mpmin; two_Mp <= two_Mpmax; two_Mp++)
			{
			  int i2 = JMarray_idx (two_Jp, two_Mp);
			  int two_k, two_kmin, two_p = two_Mp - two_M;
			  //	  double Ep;
			  gsl_complex Ap, c, eterm, mtxel;

			  GSL_SET_COMPLEX (&Ap, coef_real[i2], coef_imag[i2]);

			  if (gsl_complex_abs (Ap) < mol->coefmin)
			    continue;

			  GSL_SET_COMPLEX (&mtxel, 0.0, 0.0);

			  if ((two_J == two_Jp) && (two_p == 0))
			    two_kmin = 0;
			  else
			    two_kmin = 4;

			  for (two_k = kmin; two_k <= 4; k += 4)
			    {
			      gsl_complex Ekp = E->get (E, two_k / 2, two_p / 2);
			      double alphak0 = polarizability_get (mol->alpha, two_k / 2, 0);
			      double a = -alphak0 *
				dmtxel (two_J, 0, two_M, two_Jp, 0, two_Mp, two_k, two_p, 0);

			      gsl_complex c = gsl_complex_mul_real (Ekp, a);
			      
			      mtxel = gsl_complex_add (mtxel, c);
			    }

			  eterm = gsl_complex_polar (1.0, (EJ - EJp) * t);

			  c = gsl_complex_mul (eterm, mtxel);
			  c = gsl_complex_mul_real (c, f);
			  c = gsl_complex_mul (c, Ap);

			  /* Signs here reflect the factor of i in the TDSE. */
			  deriv_real[i1] += GSL_IMAG (c);
			  deriv_imag[i1] -= GSL_REAL (c);
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
linear_molecule_check_populations (const molecule_t * molecule, 
				   const double * coef)
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  const double *coef_r = coef;
  const double *coef_i = coef + mol->ncoef;
  int two_J;

  /* Check populations in Jmax manifold, and also Jmax - 1 manifold. */
  for (two_J = mol->two_Jmax - 2; two_J <= mol->two_Jmax; two_J += 2)
    {
      int two_M;
      double pop = 0.0;
      
      for (two_M = -two_J; two_M <= two_J; M += 2)
	{
	  int idx = JMarray_idx (J, M);
	  double re = coef_r[idx];
	  double im = coef_i[idx];
	  pop += (re * re) + (im * im); /* Population in this (J, M) state */
	}
      /* If total population in this J manifold is greater than
	 poptol, then it's considered unacceptable. */
      if (pop > mol->poptol)
	{
	  fprintf(stderr, "Population greater than poptol in 2J=%d. Population=%g, poptol=%g\n", 
		  two_J, pop, mol->poptol);
	  return -1;
	}
    }

  return 0;
}

linear_molecule_t *
linear_molecule_ctor(const double B, const int two_Jmax, const double T,
		     const double oddwt, const double evenwt, 
		     const double alpha_par, const double alpha_perp, 
		     const double coefmin, const double poptol)
/* B is rotational constant in wavenumbers.  two_Jmax is twice the
   maximum value of J to consider.

   oddwt, evenwt are the statistical weight of odd and even J levels
   respectively.

   alpha_par and alpha_perp are the parallel and perpendicular
   components of polarizability in Angstrom cubed.
   
   coefmin is the minimum value for a coefficient which we bother to
   propagate in the TDSE.

   poptol is the maximum allowed population in the Jmax and Jmax-1
   manifolds during propagation. If the population in these manifolds
   exceeds poptol the calculation will bail out.
*/
{
  double a_par, a_perp;
  linear_molecule_t * mol;
  int two_J;

  if (MEMORY_ALLOC(mol) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }
  
  /* Allocate parent structure and register dispatch functions. */
  molecule_dispatch_register((molecule_t *) mol,
			     linear_molecule_tdse_rhs,
			     linear_molecule_get_tdse_job,
			     linear_molecule_get_tdse_job_coef,
			     linear_molecule_get_tdse_job_weight,
			     linear_molecule_set_tdse_job_done,
			     linear_molecule_check_populations,
			     linear_molecule_tdse_worker_ctor,
			     linear_molecule_tdse_worker_dtor,
			     linear_molecule_get_ncoef,
			     linear_molecule_expval_ctor,
			     linear_molecule_expval_dtor,
			     linear_molecule_expval_calc,
			     linear_molecule_expval_add_weighted,
			     linear_molecule_expval_fwrite,
#ifdef BUILD_WITH_MPI
			     linear_molecule_expval_mpi_send,
			     linear_molecule_expval_mpi_recv,
			     linear_molecule_tdse_worker_mpi_send,
			     linear_molecule_tdse_worker_mpi_recv,
#endif
			     linear_molecule_dtor);
  mol->two_Jmax = two_Jmax;
  mol->poptol = poptol;
  mol->coefmin = coefmin;
  mol->oddwt = oddwt;
  mol->evenwt = evenwt;

  /* Convert to atomic units. */
  mol->B = WN_TO_AU(B);
 
  /* Caulculate kT in atomic units from T. */
  mol->kT = T_TO_KTAU (T);

  /* Calculate polarizability tensor in atomic units. */
  a_par = ANG3_TO_AU(alpha_par);
  a_perp = ANG3_TO_AU(alpha_perp);
  mol->alpha = polarizability_from_cart_ctor (a_perp, a_perp, a_par);

  if (mol->alpha == NULL)
    {
      fprintf(stderr, "Failed to allocate memory for mol->alpha.\n");
      linear_molecule_dtor((molecule_t *) mol);
      return NULL;
    }

  /* Number of complex coefficients. */
  mol->ncoef = JMarray_dim (mol->Jmax);

  /* Calculate and store the partition function. */
  mol->partfn = 0.0;
  for (two_J = 0; two_J <= mol->two_Jmax; J += 2)
    {
      double dim = two_J + 1.0;
      /* FIXME: only works for integer J. */
      double wt = GSL_IS_ODD (two_J / 2) ? mol->oddwt : mol->evenwt;
      double E = linear_molecule_energy(mol, two_J);
      mol->partfn += dim * wt * exp (-E / mol->kT);
    }

  /* Setup array for storing TDSE job status. */
  mol->job_status = JMarray_int_ctor(mol->two_Jmax);
  if (mol->job_status == NULL)
    {
      fprintf(stderr, "Failed to allocate memory for mol->tdse_status.\n");
      linear_molecule_dtor((molecule_t *)mol);
      return NULL;
    }

  for (two_J = 0; two_J <= mol->two_Jmax; two_J += 2)
    {
      int two_M;

      for (two_M = -two_J; two_M <= two_J; two_M++)
	JMarray_int_set (mol->job_status, two_J, two_M, __TODO);
    }

  return mol;
}

molecule_t *
linear_molecule_cfg_parse_ctor(const config_setting_t *cfg)
{
  linear_molecule_t * mol;
  double B, alpha_par, alpha_perp, T, coefmin, poptol;
  double oddwt, evenwt;
  int two_Jmax;

  if (!(config_setting_lookup_float (cfg, "B", &B) &&
	config_setting_lookup_float (cfg, "alpha_par", &alpha_par) &&
	config_setting_lookup_float (cfg, "alpha_perp", &alpha_perp) &&
	config_setting_lookup_float (cfg, "T", &T) &&
	config_setting_lookup_float (cfg, "oddwt", &oddwt) &&
	config_setting_lookup_float (cfg, "evenwt", &evenwt) &&
	config_setting_lookup_float (cfg, "coefmin", &coefmin) &&
	config_setting_lookup_float (cfg, "poptol", &poptol) &&
	config_setting_lookup_int (cfg, "two_Jmax", &two_Jmax) 
	))
    {
      fprintf(stderr, "Incomplete linear molecule configuration.\n");
      return NULL;
    }

  mol = linear_molecule_ctor (B, two_Jmax, T, oddwt, evenwt, alpha_par, alpha_perp, 
			      coefmin, poptol);
  if (mol == NULL)
    fprintf (stderr, "Failed to create linear_molecule structure after parsing config.\n");

  return (molecule_t *) mol;
}

void
linear_molecule_dtor(molecule_t * molecule)
{ 
  linear_molecule_t * mol = (linear_molecule_t *) molecule;

  if (mol->job_status != NULL)
    JMarray_int_dtor(mol->job_status);

  if (mol->alpha != NULL)
    polarizability_dtor(mol->alpha);

  MEMORY_FREE (mol);
}

#undef __TODO
#undef __STARTED
#undef __DONE


