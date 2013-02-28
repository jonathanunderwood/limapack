#ifndef __JKMARRAY_int_H__
#define __JKMARRAY_int_H__

typedef struct _JKMarray_int JKMarray_int_t;

JKMarray_int_t *JKMarray_int_ctor (const int Jmax);
void JKMarray_int_dtor (JKMarray_int_t * a);
int JKMarray_int_get (JKMarray_int_t * a, const int J, const int K,
		      const int M);
void JKMarray_int_set (JKMarray_int_t * a, const int J, const int K,
		       const int M, const int val);

#endif /* __JKMARRAY_int_H__ */
