#ifndef __MOLECULE_LINEAR_H__
#define __MOLECULE_LINEAR_H__

#include "molecule.h"
#include "jmarray.h"
#include "polarizability.h"

typedef struct _linear_molecule
{
  molecule_t parent;
  double B;
  double kT, partfn;
  double coefmin;
  double poptol;
  double oddwt, evenwt;
  polarizability_t * alpha;
  int Jmax;
  int ncoef;
  int nexpval;
  JMarray_int_t * job_status;
} linear_molecule_t;

linear_molecule_t * linear_molecule_ctor();
void linear_molecule_dtor(linear_molecule_t * mol);
void linear_molecule_dispatched_dtor(molecule_t * mol);
molecule_t * linear_molecule_cfg_parse_ctor(const config_setting_t *cfg);

typedef struct _linear_molecule_tdse_worker
{
  molecule_tdse_worker_t parent;
  int J, M;
} linear_molecule_tdse_worker_t;

/* These functions are registered in the molecule_t parent structure
   for dispatch/callback. */
int linear_molecule_tdse_rhs (const molecule_t *molecule, 
			      const laser_collection_t *lasers, 
			      const double t, const double *coef, 
			      double *deriv);
int linear_molecule_get_tdse_job(molecule_t *molecule, molecule_tdse_worker_t *worker,
				 double *coef,  double * weight);
void linear_molecule_set_tdse_job_done (molecule_t *molecule, 
					molecule_tdse_worker_t *worker);
int linear_molecule_check_populations (const molecule_t *molecule, 
				       const double * coef);
molecule_tdse_worker_t * linear_molecule_tdse_worker_ctor (const molecule_t *molecule);
void linear_molecule_tdse_worker_dtor (const molecule_t *molecule, 
				       molecule_tdse_worker_t * worker);
int linear_molecule_get_ncoef(const molecule_t * molecule);
int linear_molecule_get_nexpval(const molecule_t *molecule);

#endif /* __MOLECULE_LINEAR_H__ */
