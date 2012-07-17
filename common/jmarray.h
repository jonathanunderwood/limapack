#ifndef __JMARRAY_H__
#define __JMARRAY_H__

typedef struct JMarray_t
{
  double (*get) (struct JMarray_t * self, const int J, const int M);
  void (*set) (struct JMarray_t * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
} JMarray;

int JMarray_init (JMarray * a, const int Jmax);
void JMarray_free (JMarray * a);
double JMarray_get (JMarray * a, const int J, const int M);
void JMarray_set (JMarray * a, const int J, const int M, const double val);

#endif /* __JMARRAY_H__ */
