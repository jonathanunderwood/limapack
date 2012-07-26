#ifndef __JKMARRAY_H__
#define __JKMARRAY_H__

typedef struct _JKMarray
{
  double (*get) (struct _JKMarray * self, const int J, const int K,
		 const int M);
  void (*set) (struct _JKMarray * self, const int J, const int K,
	       const int M, const double val);
  double *data;
  int Jmax;
  int dim;
} JKMarray_t;

JKMarray_t * JKMarray_ctor (const int Jmax);
void JKMarray_dtor (JKMarray_t * a);
double JKMarray_get (JKMarray_t * a, const int J, const int K, const int M);
void JKMarray_set (JKMarray_t * a, const int J, const int K, const int M,
		   const double val);

#endif /* __JKMARRAY_H__ */
