#ifndef __JKMARRAY_H__
#define __JKMARRAY_H__

typedef struct JKMarray_t
{
  double (*get) (struct JKMarray_t * self, const int J, const int K,
		 const int M);
  void (*set) (struct JKMarray_t * self, const int J, const int K,
	       const int M, const double val);
  double *data;
  int Jmax;
  int dim;
} JKMarray;

int JKMarray_init (JKMarray * a, const int Jmax);
void JKMarray_free (JKMarray * a);
double JKMarray_get (JKMarray * a, const int J, const int K, const int M);
void JKMarray_set (JKMarray * a, const int J, const int K, const int M,
		   const double val);

#endif /* __JKMARRAY_H__ */
