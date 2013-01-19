#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "dcmsq.h"
#include "dcmsq_mpi.h"

/* Note that these functions assume contiguous data storage in
   dcmsq_expval_t. If this structure is ever changed to contain
   pointers, a different approach will be needed. */

int 
dcmsq_expval_mpi_send(dcmsq_expval_t *expval, 
		      int dest, int tag, MPI_Comm comm)
{
  MPI_Send(expval, sizeof(dcmsq_expval_t), MPI_BYTE, dest, tag, comm);
}

int
dcmsq_expval_mpi_recv(dcmsq_expval_t *expval, 
		      int dest, int tag, MPI_Comm comm)
{
  MPI_Send(expval, sizeof(dcmsq_expval_t), MPI_BYTE, dest, 
	   tag, comm, status);
}
