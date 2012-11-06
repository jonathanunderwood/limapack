#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include "laser.h"

typedef enum _molecule_type {
  MT_LINROT,
  MT_SYMROT,
  MT_ASYMROT
} molecule_type_t;

typedef enum _molecule_tdse_worker_state {
  TW_UNITIALIZED, TW_WAITING, TW_WORKING, TW_DEAD
} molecule_tdse_worker_state_t;

typedef enum _molecule_tdse_job_state {
  TJ_TODO, TJ_STARTED, TJ_DONE
} molecule_tdse_job_state_t;

/* Infrastructure for tracking worker process states. These will need
   to be "sub-classed" for each specific molecule type. */
typedef struct _molecule_tdse_worker
{
  molecule_tdse_worker_state_t state;
} molecule_tdse_worker_t;

molecule_tdse_worker_t *molecule_tdse_worker_ctor();
void molecule_tdse_worker_dtor(molecule_tdse_worker_t * worker);

typedef struct _molecule 
{
  /* Functions which need to be registered by the different molecule
     types. */
  int (*tdse_rhs) (const struct _molecule *mol, const laser_collection_t *lasers, 
		   const double t, const double *coef, double *deriv);
  int (*get_tdse_job) (struct _molecule *self, molecule_tdse_worker_t *worker,
		       double *coef, double *weight);
  void (*set_tdse_job_done) (struct _molecule *self, 
			     molecule_tdse_worker_t *worker);
  int (*check_populations)(const struct _molecule *self, 
			   const double * coef);
  molecule_tdse_worker_t * (*tdse_worker_ctor) (const struct _molecule *self);
  void (*tdse_worker_dtor)(const struct _molecule *self,
			   molecule_tdse_worker_t * worker);
  int (*get_ncoef) (const struct _molecule *self);
  int (*get_nexpval) (const struct _molecule *self);
  void (*dtor)(struct _molecule *self);
} molecule_t;

/* Function to calculate the right hand side of the time-dependent
   Schrodinger equation - arguments defined by the GSL ODE
   integrator. This is a generic function for all molecular types,
   which in turn dispatches tdse_rhs_mtxel which is molecule
   specific. */
int molecule_tdse_rhs (const double t, const double coefs[], 
		       double derivs[], void *params);

void
molecule_dispatch_register(molecule_t * molecule,
			   int (*tdse_rhs) (const molecule_t *mol, 
					    const laser_collection_t *lasers, 
					    const double t, const double *coef, double *deriv),
			   int (*get_tdse_job) (molecule_t *self, 
						molecule_tdse_worker_t *worker,
						double *coef, double *weight),
			   void (*set_tdse_job_done) (molecule_t *self, 
						      molecule_tdse_worker_t *worker),
			   int (*check_populations)(const molecule_t *self, 
						    const double * coef),
			   molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
			   void (*tdse_worker_dtor)(const molecule_t *self, 
						    molecule_tdse_worker_t * worker),
			   int (*get_ncoef) (const molecule_t *self),
			   int (*get_nexpval) (const molecule_t *self),
			   void (*dtor)(molecule_t *self)
			   );

molecule_t *
molecule_ctor(int (*tdse_rhs) (const molecule_t *mol, const laser_collection_t *lasers, 
			       const double t, const double *coef, double *deriv),
	      int (*get_tdse_job) (molecule_t *self, 
				   molecule_tdse_worker_t *worker,
				   double *coef, double *weight),
	      void (*set_tdse_job_done) (molecule_t *self, 
					 molecule_tdse_worker_t *worker),
	      int (*check_populations)(const molecule_t *self, 
				       const double * coef),
	      molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
	      void (*tdse_worker_dtor)(const molecule_t *self, 
				       molecule_tdse_worker_t * worker),
	      int (*get_ncoef) (const molecule_t *self),
	      int (*get_nexpval) (const molecule_t *self),
	      void (*dtor)(struct _molecule *self)
	      );
void molecule_dtor(molecule_t * molecule);

/* Generic configuration parser - this calls the molecule type
   specific parsers as required. */
molecule_t * molecule_cfg_parse_ctor (const config_t * cfg);

#endif /* __MOLECULE_H__ */
