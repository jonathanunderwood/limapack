/* Matrix elements of the squared direction cosine matrix elements */
#ifndef __DCMSQ_H__
#define __DCMSQ_H__

#include <gsl/gsl_complex.h>
#include <hdf5.h>

typedef enum dcmsq_axes
{ Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz } 
dcmsq_axes_t;

typedef struct _dcmsq_mtxel dcmsq_mtxel_t;
typedef struct _dcmsq_expval dcmsq_expval_t;

dcmsq_mtxel_t * dcmsq_mtxel_ctor ();
void dcmsq_mtxel_dtor (dcmsq_mtxel_t * dcmsq_mtxel);

dcmsq_expval_t * dcmsq_expval_ctor ();
void dcmsq_expval_dtor (dcmsq_expval_t * dcmsq_expval);

void dcmsq_mtxel_calc (const int J, const int K, const int M,
		       const int Jp, const int Kp, const int Mp, 
		       dcmsq_mtxel_t *mtxel);

void dcmsq_expval_add_mtxel_weighted_complex (dcmsq_expval_t *expval, 
					      const dcmsq_mtxel_t *mtxel,
					      const gsl_complex weight);

void dcmsq_expval_add_weighted (dcmsq_expval_t *expval_a, 
				const dcmsq_expval_t *expval_b,
				const double weight);

void dcmsq_expval_zero(dcmsq_expval_t *dcmsq);

int dcmsq_fwrite(const dcmsq_expval_t *expval, const hid_t *location);

#ifdef BUILD_WITH_MPI
#include <mpi.h>
int dcmsq_expval_mpi_send(const dcmsq_expval_t *expval, 
			  int dest, int tag, MPI_Comm comm);

int dcmsq_expval_mpi_recv(dcmsq_expval_t *expval, 
			  int dest, int tag, MPI_Comm comm);
#endif

//double dcmsq_get(const dcmsq_buff_t *buff, const dcmsq_axes_t axes);

double dcmsq_Xx (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Xy (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Xz (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Yx (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Yy (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Yz (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Zx (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Zy (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);

double dcmsq_Zz (const int J, const int K, const int M,
		 const int Jp, const int Kp, const int Mp);
#endif /* __DCMSQ_H__ */
