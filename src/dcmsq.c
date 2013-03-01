#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <hdf5.h>

#include "dmtxel.h"
#include "dcmsq.h"
#include "hdf5_gsl_complex.h"
#include "memory.h"

struct _dcmsq_mtxel
{
  double data[9];
};

struct _dcmsq_expval
{
  gsl_complex data[9];
};

dcmsq_expval_t *
dcmsq_expval_ctor ()
{
  dcmsq_expval_t * dcmsq = NULL;

  if (MEMORY_ALLOC (dcmsq) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  return dcmsq;
}

void
dcmsq_expval_dtor (dcmsq_expval_t * dcmsq)
{
  MEMORY_FREE (dcmsq);
}

dcmsq_mtxel_t * 
dcmsq_mtxel_ctor ()
{
  dcmsq_mtxel_t * dcmsq = NULL;

  if (MEMORY_ALLOC (dcmsq) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  return dcmsq;
}

void 
dcmsq_mtxel_dtor (dcmsq_mtxel_t * dcmsq)
{
  MEMORY_FREE (dcmsq);
}

void
dcmsq_mtxel_calc (dcmsq_mtxel_t * mtxel, 
		  const int J, const int K, const int M,
		  const int Jp, const int Kp, const int Mp)
/* Calculates the matrix elements of the squared direction cosine
   matrix elements <J K M| <A|b>^2 |Jp Kp Mp> for A = X, Y, Z and b =
   x, y, z. Integer angular momentum only at this point. */
{
  dcmsq_axes_t idx;
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double f = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 2);
  double g = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, -2);
  double h = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 2);
  double i = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  /* <J K Mx| <X|x>^2 |Jp Kp Mp> */
  mtxel->data[idx = Xx] = (j / 3.0) + a / 6.0 + (f + g + h + i) / 4.0 -
    (b + c + d + e) / (2.0 * sqrt (6.0));

  /* <J K M| <X|y>^2 |Jp Kp Mp> */
  mtxel->data[idx = Xy] = (j / 3.0) + a / 6.0 - (f + g + h + i) / 4.0 -
    (d + e - b - c) / (2.0 * sqrt (6.0));

  /* <J K M| <X|z>^2 |Jp Kp Mp> */
  mtxel->data[idx = Xz] = (j - a) / 3.0 + (d + e) / sqrt (6.0);

  /* <J K M| <Y|x>^2 |Jp Kp Mp> */
  mtxel->data[idx = Yx] = (j / 3.0) + a / 6.0 - (f + g + h + i) / 4.0 -
    (b + c - d - e) / (2.0 * sqrt (6.0));

  /* <J K M| <Y|y>^2 |Jp Kp Mp> */
  mtxel->data[idx = Yy] = (j / 3.0) + a / 6.0 + (f + g + h + i) / 4.0 +
    (b + c - d - e) / (2.0 * sqrt (6.0));

  /* <J K M| <Y|z>^2 |Jp Kp Mp> */
  mtxel->data[idx = Yz] = (j - a) / 3.0 - (d + e) / sqrt (6.0);

  /* <J K M| <Z|x>^2 |Jp Kp Mp> */
  mtxel->data[idx = Zx] = (j - a) / 3.0 + (c - b) / sqrt (6.0);

  /* <J K M| <Z|y>^2 |Jp Kp Mp> */
  mtxel->data[idx = Zy] = (j - a) / 3.0 - (b + c) / sqrt (6.0);

  /* <J K M| <Z|z>^2 |Jp Kp Mp> */
  mtxel->data[idx = Zz] = (j + 2.0 * a) / 3.0;

  return;
}

void dcmsq_expval_add_mtxel_weighted_complex (dcmsq_expval_t *expval, 
					      const dcmsq_mtxel_t *mtxel,
					      const gsl_complex weight)
/* For each expectation value add to it the matrix element in mtxel
   multiplied by the complex coefficient weight. 

   i.e expval = expval + weight * mtxel. */
{
  int i;

  for (i = 0; i < 9; i++)
    {
      gsl_complex weighted_mtxel = 
	gsl_complex_mul_real (weight, mtxel->data[i]);

      expval->data[i] = 
	gsl_complex_add (expval->data[i], weighted_mtxel);
    }
  return;
}

void dcmsq_expval_add_weighted (dcmsq_expval_t *a, 
				const dcmsq_expval_t *b,
				const double weight)
/* Add expectation values as a = a + weight * b */
{
  int i;

  for (i = 0; i < 9; i++)
    {
      gsl_complex c = gsl_complex_mul_real (b->data[i], weight);
      a->data[i] = gsl_complex_add (a->data[i], c);
    }

  return;
}

void 
dcmsq_expval_zero(dcmsq_expval_t *dcmsq)
{
  int i;
  gsl_complex zero = gsl_complex_rect (0.0, 0.0);

  for (i = 0; i < 9; i++)
    dcmsq->data[i] = zero;
}

#ifdef BUILD_WITH_MPI
#include <mpi.h>
/* Note that these functions assume contiguous data storage in
   dcmsq_expval_t. If this structure is ever changed to contain
   pointers, a different approach will be needed. */

int
dcmsq_expval_mpi_send(const dcmsq_expval_t *expval, 
		      int dest, int tag, MPI_Comm comm)
{
  /* Note that the MPI_Send function (and the whole of the MPI 2 API)
     doesn't use the const qualifier, so we have to recast to void
     here to drop the const. This is idiotic, but there ya go. This
     is fixed in the MPI 3 spec. */
  return MPI_Send((void *) expval, sizeof(dcmsq_expval_t), MPI_BYTE, dest, tag, comm);
}

int
dcmsq_expval_mpi_recv(dcmsq_expval_t *expval, 
		      int dest, int tag, MPI_Comm comm)
{
  return MPI_Recv(expval, sizeof(dcmsq_expval_t), MPI_BYTE, dest, 
		  tag, comm, MPI_STATUS_IGNORE);
}
#endif

int
dcmsq_fwrite(const dcmsq_expval_t *expval, const hid_t *location)
{
  hid_t cmplx_type;
  int status;
  /* Data space variables. */
  hid_t dataspace_id;
  int rank = 1;
  hsize_t dim[1] = {9}; 
  /* Data set variables. */
  hid_t dataset_id;

  /* Define and hdf5 compound data type to represent a gsl_complex type. */
  cmplx_type = hdf5_gsl_complex_type_create();
  if (cmplx_type < 0)
    {
      fprintf(stderr, "%s %d: Failed to create HDF5 complex type.\n", __func__, __LINE__);
      return -1;
    }

  /* Create dataspace. */
  dataspace_id = H5Screate_simple(rank, dim, NULL);
  if (dataspace_id < 0)
    {
      fprintf(stderr, "%s %d: Failed to create dataspace.\n", __func__, __LINE__);
      return -1;
    }

  /* Create dataset. Note that we don't have a leading "/" as this
     would force the data set into the root of the file, whereas
     location may refer to a group. */
  dataset_id = H5Dcreate(*location, "dcmsq", cmplx_type, dataspace_id, 
			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0)
    {
      fprintf(stderr, "%s %d: Failed to create dataset.\n", __func__, __LINE__);
      return -1;
    }

  /* Write to dataset. */
  status = H5Dwrite (dataset_id, cmplx_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, expval);
  if (status)
    {
      fprintf(stderr, "%s %d: Failed to write to dataset.\n", __func__, __LINE__);
      return -1;
    }
  /* Cleanup. */
  status += H5Dclose (dataset_id);
  status += H5Sclose (dataspace_id);
  status += H5Tclose (cmplx_type);

  return status;
}

double
dcmsq_Xx (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <X|x>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double f = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 2);
  double g = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, -2);
  double h = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 2);
  double i = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j / 3.0) + a / 6.0 + (f + g + h + i) / 4.0 -
    (b + c + d + e) / (2.0 * sqrt (6.0));
}

double
dcmsq_Xy (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <X|y>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double f = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 2);
  double g = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, -2);
  double h = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 2);
  double i = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j / 3.0) + a / 6.0 - (f + g + h + i) / 4.0 -
    (d + e - b - c) / (2.0 * sqrt (6.0));
}

double
dcmsq_Xz (const int J, const int K, const int M, const int Jp, const int Kp,
	  const int Mp)
     /* <J K M| <X|z>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j - a) / 3.0 + (d + e) / sqrt (6.0);
}

double
dcmsq_Yx (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <Y|x>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double f = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 2);
  double g = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, -2);
  double h = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 2);
  double i = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j / 3.0) + a / 6.0 - (f + g + h + i) / 4.0 -
    (b + c - d - e) / (2.0 * sqrt (6.0));

}

double
dcmsq_Yy (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <Y|y>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double f = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 2);
  double g = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, -2);
  double h = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 2);
  double i = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j / 3.0) + a / 6.0 + (f + g + h + i) / 4.0 +
    (b + c - d - e) / (2.0 * sqrt (6.0));
}

double
dcmsq_Yz (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
	  /* <J K M| <Y|z>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double d = dmtxel (J, K, M, Jp, Kp, Mp, 2, 2, 0);
  double e = dmtxel (J, K, M, Jp, Kp, Mp, 2, -2, 0);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j - a) / 3.0 - (d + e) / sqrt (6.0);
}

double
dcmsq_Zx (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <Z|x>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j - a) / 3.0 + (c - b) / sqrt (6.0);
}

double
dcmsq_Zy (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <Z|y>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double b = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 2);
  double c = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, -2);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j - a) / 3.0 - (b + c) / sqrt (6.0);
}

double
dcmsq_Zz (const int J, const int K, const int M,
	  const int Jp, const int Kp, const int Mp)
     /* <J K M| <Z|z>^2 |Jp Kp Mp> */
{
  double a = dmtxel (J, K, M, Jp, Kp, Mp, 2, 0, 0);
  double j = dmtxel (J, K, M, Jp, Kp, Mp, 0, 0, 0);

  return (j + 2.0 * a) / 3.0;
}
