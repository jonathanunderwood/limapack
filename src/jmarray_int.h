#ifndef __JMARRAY_int_H__
#define __JMARRAY_int_H__

typedef struct _JMarray_int
{
  double (*get) (struct _JMarray_int * self, const int J, const int M);
  void (*set) (struct _JMarray_int * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
} JMarray_int_t;

JMarray_int_t * JMarray_int_ctor (const int Jmax);
void JMarray_int_dtor (JMarray_int_t * a);
int JMarray_int_get (JMarray_int_t * a, const int J, const int M);
void JMarray_int_set (JMarray_int_t * a, const int J, const int M, const int val);

#endif /* __JMARRAY_int_H__ */
