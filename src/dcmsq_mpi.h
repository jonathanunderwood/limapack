#ifndef __DCMSQ_MPI_H__
#define __DCMSQ_MPI_H__

#include <mpi.h>
#include "dcmsq.h"

int dcmsq_expval_mpi_send(dcmsq_expval_t *expval, 
			  int dest, int tag, MPI_Comm comm );

int
dcmsq_expval_mpi_recv(dcmsq_expval_t *expval, 
		      int dest, int tag, MPI_Comm comm );

#endif
