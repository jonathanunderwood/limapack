#ifndef __DMTXEL_H__
#define __DMTXEL_H__

/* Matrix elements of the Wigner rotation matrix elements. Arguments
   passed are twice the angular momenta. */
double dmtxel (const int two_J, const int two_K, const int two_M,
	       const int two_Jp, const int two_Kp, const int two_Mp,
	       const int two_k, const int two_p, const int two_q);

#endif /* __DMTXEL_H__ */
