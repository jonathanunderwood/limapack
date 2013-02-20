#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include <hdf5.h>

#ifdef BUILD_WITH_MPI
#include <mpi.h>
#endif

#include "laser.h"

typedef enum _molecule_type {
  MT_LINROT,
  MT_SYMROT,
  MT_ASYMROT
} molecule_type_t;

/* This struct will keep track of what specific TDSE propagation job a
   process is working on. It will need to be "sub-classed" by each
   molecular model which will implement a struct with a
   molecule_tdse_worker_t as its first member:

   typedef struct _my_molecule_tdse_worker
   {
     molecule_tdse_worker_t parent;
     int J, K, M, a;
   } my_molecule_tdse_worker_t;

*/
#define MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH 128
typedef struct _molecule_tdse_worker
{
  char description[MOLECULE_TDSE_WORKER_DESCRIPTION_LENGTH];
} molecule_tdse_worker_t;

/* The following typedefs are actually incomplete - nowhere do we
   define the corresponding structures. Rather, for each molecular
   model we define specific structs for that model and recast as
   necessary. We elect to use incomplete definitions rather than void
   here to give us some level of type checking. As such they are
   "placeholder" structs.  The first represents a struct that is used
   to store expectation values for each time step. */ 
typedef struct _molecule_expval molecule_expval_t;

/* The molecule_t struct is simply a long list of call back function
   pointers that each molecular model needs to register. Each
   molecular model should implement its own specific structure with a
   molecule_t struct as its first member:

   struct my_new_molecule_model_t
   {
     molecule_t parent;
     int some_param;
     double some_other_param;
     ...
   };
   
   As such, a my_new_molecule_model_t pointer can be recast as a
   molecule_t pointer at any point. In this way we can write general
   functions for all molecular models as long as they provide
   conformant callback functions in their struct for later
   dispatch. See eg. odesys.c. We provide a function here for
   registering the callback functions as well,, and so each model
   should make use of this. */
  
typedef struct _molecule 
{
#ifdef BUILD_WITH_PTHREADS
  pthread_mutex_t lock;
#endif
  int (*tdse_rhs) (const struct _molecule *self, const laser_collection_t *lasers, 
		   const double t, const double *coef, double *deriv);
  int (*get_tdse_job) (struct _molecule *self, molecule_tdse_worker_t *worker);
  void (*get_tdse_job_coef) (const struct _molecule *self, 
			     const molecule_tdse_worker_t *worker,
			     double *coef);
  double (*get_tdse_job_weight) (const struct _molecule *self, 
				 const molecule_tdse_worker_t *worker);
  void (*set_tdse_job_done) (struct _molecule *self, 
			     molecule_tdse_worker_t *worker);
  int (*check_populations)(const struct _molecule *self, 
			   const double * coef);
  molecule_tdse_worker_t * (*tdse_worker_ctor) (const struct _molecule *self);
  void (*tdse_worker_dtor)(const struct _molecule *self,
			   molecule_tdse_worker_t * worker);
  int (*get_ncoef) (const struct _molecule *self);
  molecule_expval_t * (*expval_ctor) (const struct _molecule *self);
  void (*expval_dtor) (const struct _molecule *self, molecule_expval_t *expval);
  int (*expval_calc)  (const struct _molecule *self, const double *coef, 
			const double t, molecule_expval_t *expval);
  void (*expval_add_weighted) (const struct _molecule *self, 
			       molecule_expval_t *a, const molecule_expval_t *b,
			       const double weight);
  int (*expval_fwrite) (const struct _molecule *self,
			const molecule_expval_t *expval,
			const hid_t *location);
#ifdef BUILD_WITH_MPI
  int (*expval_mpi_send) (const struct _molecule *self, 
			  const molecule_expval_t *molecule, 
			  int dest, int tag, MPI_Comm comm);
  int (*expval_mpi_recv) (const struct _molecule *self, 
			  molecule_expval_t *molecule, 
			  int dest, int tag, MPI_Comm comm);
  int (*tdse_worker_mpi_send) (const struct _molecule *self, 
			       const molecule_tdse_worker_t *worker,
			       int dest, int tag, MPI_Comm comm);
  int (*tdse_worker_mpi_recv) (const struct _molecule *self, 
			       molecule_tdse_worker_t *worker,
			       int dest, int tag, MPI_Comm comm);
#endif
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
						molecule_tdse_worker_t *worker),
			   void (*get_tdse_job_coef) (const molecule_t *self, 
						      const molecule_tdse_worker_t *worker,
						      double *coef),
			   double (*get_tdse_job_weight) (const molecule_t *self, 
							  const molecule_tdse_worker_t *worker),
			   void (*set_tdse_job_done) (molecule_t *self, 
						      molecule_tdse_worker_t *worker),
			   int (*check_populations)(const molecule_t *self, 
						    const double * coef),
			   molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
			   void (*tdse_worker_dtor)(const molecule_t *self, 
						    molecule_tdse_worker_t * worker),
			   int (*get_ncoef) (const molecule_t *self),
			   molecule_expval_t * (*expval_ctor) (const struct _molecule *self),
			   void (*expval_dtor) (const struct _molecule *self, 
						molecule_expval_t *expval),
			   int (*expval_calc)  (const molecule_t *self, const double *coef,
						 const double t, molecule_expval_t *expval),
			   void (*expval_add_weighted) (const molecule_t *self, 
							molecule_expval_t *a, 
							const molecule_expval_t *b,
							const double weight),
			   int (*expval_fwrite) (const molecule_t *molecule,
						 const molecule_expval_t *expvalue,
						 const hid_t *location),
#ifdef BUILD_WITH_MPI
			   int (*expval_mpi_send) (const molecule_t *self, 
						   const molecule_expval_t *expval, 
						   int dest, int tag, MPI_Comm comm),
			   int (*expval_mpi_recv) (const molecule_t *self, 
						   molecule_expval_t *molecule, 
						   int dest, int tag, MPI_Comm comm),
			   int (*tdse_worker_mpi_send) (const molecule_t *self, 
							const molecule_tdse_worker_t *worker,
							int dest, int tag, MPI_Comm comm),
			   int (*tdse_worker_mpi_recv) (const molecule_t *self, 
							molecule_tdse_worker_t *worker,
							int dest, int tag, MPI_Comm comm),
#endif
			   void (*dtor)(molecule_t *self)
			   );

/* Generic configuration parser - this calls the molecule type
   specific parsers as required. */
molecule_t * molecule_cfg_parse_ctor (const config_t * cfg);

#endif /* __MOLECULE_H__ */
