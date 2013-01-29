#ifndef __ODESYS_H__
#define __ODESYS_H__

/* The intention of this module is to hide as much as possible the
   details of the ODE solver - in the future may use something other
   than GSL, but the API should remain the more or less the same. */

#include <gsl/gsl_odeiv.h>

/* The integration function should return ODESYS_SUCCESS. */
#include <gsl/gsl_errno.h>
#define ODESYS_SUCCESS GSL_SUCCESS

typedef struct _odesys odesys_t;

odesys_t * odesys_parse_cfg_from_file_ctor (const char *file, 
					    const int max_cfg_file_size);
odesys_t * odesys_parse_cfg_from_buffer_ctor (const char *buffer);

void odesys_dtor (odesys_t * ode);
void odesys_reset (odesys_t * ode);
int odesys_step (odesys_t * ode, 
		 const double t1, const double t2, 
		 double *coef);

int odesys_tdse_propagate_simple (odesys_t *odesys);

int odesys_expval_fwrite(const odesys_t *ode, const char *filename);

#ifdef BUILD_WITH_MPI
int odesys_tdse_propagate_mpi_master (odesys_t *odesys);
int odesys_tdse_propagate_mpi_slave (odesys_t *odesys);
#endif /* BUILD_WITH_MPI */

#endif /* __ODESYS_H__ */

