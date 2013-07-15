#ifndef __JKMARRAY_@@TYPE@@_H__
#define __JKMARRAY_@@TYPE@@_H__

typedef struct _JKMarray_@@TYPE@@ JKMarray_@@TYPE@@_t;

JKMarray_@@TYPE@@_t * JKMarray_@@TYPE@@_ctor (const int two_Jmax);
void JKMarray_@@TYPE@@_dtor (JKMarray_@@TYPE@@_t * a);
@@TYPE@@ JKMarray_@@TYPE@@_get (JKMarray_@@TYPE@@_t * a, const int two_J, const int two_K, 
				const int two_M);
void JKMarray_@@TYPE@@_set (JKMarray_@@TYPE@@_t * a, const int two_J, const int two_K, 
			    const int two_M, const @@TYPE@@ val);

#endif /* __JKMARRAY_@@TYPE@@_H__ */
