#ifndef __AU_H__
#define __AU_H__

/* Some macros to convert to/from atomic units. */

/* Convert from wavenumbers to atomic units of energy. */
#define WN_TO_AU(A) ((A) * 4.556334874e-6)

/* Convert from ps to atomic units of time. */
#define PS_TO_AU(A) ((A) * 4.134137334e4)

/* Convert from atomic units of time to ps. */
#define AU_TO_PS(A) ((A) * 2.418884327e-5)

/* Convert T in K to the value of kT in atomic units */
#define T_TO_KTAU(T) ((T) * 3.166832608e-6)

/* Convert polarizability in Angstrom cubed to atomic units. */
#define ANG3_TO_AU(A) ((A) *6.74833499)

/* Convert electric field in V/m to atomic units. */
#define E_TO_AU(E) ((E) * 1.944690505e-12)

/* Convert frequency in Hertz to atomic units. */
#define HZ_TO_AU(HZ) ((HZ) / 4.13413732e16)


#endif /* __AU_H__ */
