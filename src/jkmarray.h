#ifndef __JKMARRAY_H__
#define __JKMARRAY_H__

/* These two utility functions are defined here to allow for use
   outside of the JKMarray_foo_t types. The first maps a (J, K, M)
   triplet, with K=-J..J and M=-J..J, to an index in a 1D array. The
   second gives the dimension needed for that array given a maximum
   value of J. */
extern inline int
JKMarray_idx (const int J, const int K, const int M)
{
  return (((4 * J * J + 5) * J) / 3 + 2 * J * J + K * (2 * J + 1) + M);
}

extern inline int
JKMarray_dim (const int Jmax)
{
  return (((4 * Jmax * Jmax + 11) * Jmax) / 3 + (4 * Jmax * Jmax) + 1);
}

#endif /* __JKMARRAY_H__ */
