#ifndef __JMARRAY_H__
#define __JMARRAY_H__

typedef struct _JMarray
{
  double (*get) (struct _JMarray * self, const int J, const int M);
  void (*set) (struct _JMarray * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
} JMarray_t;

JMarray_t * JMarray_ctor (const int Jmax);
void JMarray_dtor (JMarray_t * a);
double JMarray_get (JMarray_t * a, const int J, const int M);
void JMarray_set (JMarray_t * a, const int J, const int M, const double val);

#endif /* __JMARRAY_H__ */
