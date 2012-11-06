#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>

#include <libconfig.h>

#include "molecule.h"
#include "laser.h"
#include "polarizability.h"
#include "jmarray.h"
#include "memory.h"
#include "molecule_linear.h"
#include "dmtxel.h"
#include "au.h"

#define __LM_NEXPVAL__ 9


static inline double
linear_molecule_energy(const linear_molecule_t *mol, const int J)
{
  return mol->B * J * (J + 1);
}

static double
linear_molecule_boltzmann_statwt(const linear_molecule_t * mol,
				 const int J)
/* Returns the statistical weight of a single (J,M) state in a
   Boltzmann distribution, hence no (2J+1) degeneracy pre-factor. */
{ 
  double wt = GSL_IS_ODD (J) ? mol->oddwt : mol->evenwt; 
  double E = linear_molecule_energy (mol, J); 

  return wt * exp (-E / mol->kT) / mol->partfn; 
}

int
linear_molecule_get_nexpval(const molecule_t *molecule)
{
  return __LM_NEXPVAL__;
}

int
linear_molecule_get_tdse_job(molecule_t *molecule, molecule_tdse_worker_t *worker,
			     double *coef,  double * weight)
{
  linear_molecule_t *mol = (linear_molecule_t *) molecule;
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;
  int J, M;

  for (J = 0; J <= mol->Jmax; J++)
    for (M = -J; M <= J; M++)
      {
	molecule_tdse_job_state_t status = 
	  mol->job_status->get(mol->job_status, J, M);
	if (status == TJ_TODO)
	  {
	    double wt = linear_molecule_boltzmann_statwt(mol, J);
	    fprintf(stdout, "J: %d M: %d wt: %g\n", J, M, wt); 
	    if (wt < mol->poptol) /* Population in this state negligible */
		mol->job_status->set(mol->job_status, J, M, TJ_DONE);
	    else
	      {
		int i;

		mol->job_status->set(mol->job_status, J, M, TJ_STARTED);
		*weight = wt;

		/* Set all real and imaginary parts of the
		   wavefunction coefficients to zero. */ 
		for (i = 0; i < 2 * molecule->get_ncoef(molecule); i++)
		  coef[i] = 0.0;
		
		/* And now set the real part of the (J, M) coefficient to 1. */
		coef[JMarray_idx (J, M)] = 1.0; /* Im part is 0.0, set above. */

		w->J = J;
		w->M = M;

		return 0;
	      }
	  }
      }

  /* No jobs in the TODO state. */
  return -1;
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

  worker->parent.state = TW_WAITING;

  /* worker->parent = molecule_tdse_worker_ctor(); */
  /* if (worker->parent == NULL) */
  /*   { */
  /*     MEMORY_FREE(worker); */
  /*     return NULL; */
  /*   } */

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

void
linear_molecule_set_tdse_job_done (molecule_t *molecule, 
				   molecule_tdse_worker_t *worker)
{
  linear_molecule_t *m = 
    (linear_molecule_t *) molecule;
  linear_molecule_tdse_worker_t *w = 
    (linear_molecule_tdse_worker_t *) worker;

  m->job_status->set(m->job_status, w->J, w->M, TJ_DONE);
  w->parent.state = TW_WAITING;
}

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
	  int J;//, idx1 = -1;

	  if (env2 < 0.0)
	    continue;

	  f = env1 * env2;

	  laser_polzn_tensor_set_from_vectors
	    (E, lasers->laser[ilas]->get_polzn_vector(lasers->laser[ilas], t),
	     lasers->laser[jlas]->get_polzn_vector(lasers->laser[jlas], t));

	  for (J = 0; J <= mol->Jmax; J++)
	    {
	      int Jpmin = J - 2, Jpmax = J + 2, M;
	      double EJ;

	      if (Jpmin < 0)
		Jpmin = 0;

	      if (Jpmax > mol->Jmax)
		Jpmax = mol->Jmax;

	      EJ = linear_molecule_energy (mol, J);
	      
	      for (M = -J; M <= J; M++)
		{
		  int Jp;
		  int i1 = JMarray_idx (J, M);
		  
		  for (Jp = Jpmin; Jp <= Jpmax; Jp++)
		    {
		      int Mp, Mpmin = GSL_MAX (-Jp, M - 2);
		      int Mpmax = GSL_MIN (Jp, M + 2);
		      double EJp = linear_molecule_energy (mol, Jp);

		      for (Mp = Mpmin; Mp <= Mpmax; Mp++)
			{
			  int i2 = JMarray_idx (Jp, Mp);
			  int k, kmin, p = Mp - M;
			  //	  double Ep;
			  gsl_complex Ap, c, eterm, mtxel;

			  GSL_SET_COMPLEX (&Ap, coef_real[i2], coef_imag[i2]);

			  if (gsl_complex_abs (Ap) < mol->coefmin)
			    continue;

			  GSL_SET_COMPLEX (&mtxel, 0.0, 0.0);

			  if ((J == Jp) && (p == 0))
			    kmin = 0;
			  else
			    kmin = 2;

			  for (k = kmin; k <= 2; k += 2)
			    {
			      gsl_complex Ekp = E->get (E, k, p);
			      double a = -mol->alpha->get (mol->alpha, k, 0) *
				dmtxel (J, 0, M, Jp, 0, Mp, k, p, 0);

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
  int J;

  /* Check populations in Jmax manifold, and also Jmax - 1 manifold. */
  for (J = mol->Jmax - 1; J <= mol->Jmax; J++)
    {
      int M;
      double re, im, pop = 0.0;
      
      for (M = -J; M <= J; M++)
	{
	  int idx = JMarray_idx (J, M);
	  re = coef_r[idx];
	  im = coef_i[idx];
	  pop += (re * re) + (im * im); /* Population in this (J, M) state */
	}
      /* If total population in this J manifold is greater than
	 poptol, then it's considered unacceptable. */
      if (pop > mol->poptol)
	{
	  fprintf(stderr, "Population greater than poptol in J=%d. Population=%g, poptol=%g\n", 
		  J, pop, mol->poptol);
	  return -1;
	}
    }

  return 0;
}

linear_molecule_t *
linear_molecule_ctor(const double B, const int Jmax, const double T,
		     const double oddwt, const double evenwt, 
		     const double alpha_par, const double alpha_perp, 
		     const double coefmin, const double poptol)
/* B is rotational constant in wavenumbers.  Jmax is the maximum value
   of J to consider.  

   oddwt, evenwt are the statistical weight of odd and even J levels
   respectively.

   alpha_par and alpha_perp are the parallel and perpendicular
   components of polarizability in Angstrom cubed.
   
   coefmin is the minimum value for a coefficient which we bother to
   propagate in the TDSE.

   poptol is the minimum initial population we bother to propagate in the TDSE.
*/
{
  double a_par, a_perp;
  linear_molecule_t * mol;
  int J;

  if (MEMORY_ALLOC(mol) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }
  
  /* Allocate parent structure and register dispatch functions. */
  molecule_dispatch_register((molecule_t *) mol,
			     linear_molecule_tdse_rhs,
			     linear_molecule_get_tdse_job,
			     linear_molecule_set_tdse_job_done,
			     linear_molecule_check_populations,
			     linear_molecule_tdse_worker_ctor,
			     linear_molecule_tdse_worker_dtor,
			     linear_molecule_get_ncoef,
			     linear_molecule_get_nexpval,
			     linear_molecule_dispatched_dtor);
  mol->Jmax = Jmax;
  mol->poptol = poptol;
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
      linear_molecule_dtor(mol);
      return NULL;
    }

  /* Number of complex coefficients. */
  mol->ncoef = JMarray_dim (mol->Jmax);

  /* Calculate and store the partition function. */
  mol->partfn = 0.0;
  for (J = 0; J <= mol->Jmax; J++)
    {
      double dim = 2.0 * J + 1.0;
      double wt = GSL_IS_ODD (J) ? mol->oddwt : mol->evenwt;
      double E = linear_molecule_energy(mol, J);
      mol->partfn += dim * wt * exp (-E / mol->kT);
    }

  /* Setup array for storing TDSE job status. */
  mol->job_status = JMarray_int_ctor(mol->Jmax);
  if (mol->job_status == NULL)
    {
      fprintf(stderr, "Failed to allocate memory for mol->tdse_status.\n");
      linear_molecule_dtor(mol);
      return NULL;
    }

  return mol;
}

molecule_t *
linear_molecule_cfg_parse_ctor(const config_setting_t *cfg)
{
  linear_molecule_t * mol;
  double B, alpha_par, alpha_perp, T, coefmin, poptol;
  double oddwt, evenwt;
  int Jmax;

  if (!(config_setting_lookup_float (cfg, "B", &B) &&
	config_setting_lookup_float (cfg, "alpha_par", &alpha_par) &&
	config_setting_lookup_float (cfg, "alpha_perp", &alpha_perp) &&
	config_setting_lookup_float (cfg, "T", &T) &&
	config_setting_lookup_float (cfg, "oddwt", &oddwt) &&
	config_setting_lookup_float (cfg, "evenwt", &evenwt) &&
	config_setting_lookup_float (cfg, "coefmin", &coefmin) &&
	config_setting_lookup_float (cfg, "poptol", &poptol) &&
	config_setting_lookup_int (cfg, "Jmax", &Jmax) 
	))
    {
      fprintf(stderr, "Incomplete linear molecule configuration.\n");
      return NULL;
    }

  mol = linear_molecule_ctor (B, Jmax, T, oddwt, evenwt, alpha_par, alpha_perp, 
			      coefmin, poptol);
  if (mol == NULL)
    fprintf (stderr, "Failed to create linear_molecule structure after parsing config.\n");

  return (molecule_t *) mol;
}

void
linear_molecule_dispatched_dtor(molecule_t * mol)
/* This is a destructor function that aceepts a generic molecule_t,
   recasts and calls the true destructor - this is needed for dispatch. */
{ 
linear_molecule_dtor ((linear_molecule_t *) mol); 
}

void
linear_molecule_dtor(linear_molecule_t * mol)
{
  if (mol->job_status != NULL)
    JMarray_int_dtor(mol->job_status);

  if (mol->alpha != NULL)
    polarizability_dtor(mol->alpha);

  MEMORY_FREE (mol);
}

#undef __LM_NEXPVAL__
