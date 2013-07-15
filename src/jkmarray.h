#ifndef __JKMARRAY_H__
#define __JKMARRAY_H__

/* These utility functions are defined here to allow for use outside
   of the JKMarray_foo_t types. The first maps a (J, K, M) triplet,
   with K=-J..J and M=-J..J, to an index in a 1D array. In this array
   mapping the last index (M) is assumed to vary fastest, then the
   second (K), then the first (J). The second frunction gives the
   dimension needed for that array given a maximum value of
   J. The same functions apply for both integer and half integer J. */
extern inline int
JKMarray_idx (const int two_J, const int two_K, const int two_M)
{
  return ((two_J + 1) * (two_J  * (two_J + 2) + 3 * two_K) + 
	  3 * (two_J + two_M)) / 6;
}

extern inline int
JKMarray_dim (const int two_Jmax)
{
  return JKMarray_idx (two_Jmax, two_Jmax, two_Jmax) + 1;
}

#endif /* __JKMARRAY_H__ */
