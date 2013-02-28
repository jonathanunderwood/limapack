#ifndef __MOLECULE_ASYMROT_EIGSYS_H__
#define __MOLECULE_ASYMROT_EIGSYS_H__

typedef struct _asymrot_eigsys asymrot_eigsys_t;

asymrot_eigsys_t * asymrot_eigsys_ctor ();
void asymrot_eigsys_dtor ();

asymrot_eigsys_t *asymrot_eigsys_ctor (const int Jmax, const double Bx, 
					 const double By, const double Bz);
void asymrot_eigsys_dtor (asymrot_eigsys_t *eigsys);
double asymrot_eigsys_eigvec_coef_get (asymrot_eigsys_t *eigsys, const int J, 
				       const int n, const int K);
double asymrot_eigsys_eigval_get (asymrot_eigsys_t *eigsys, const int J, 
				  const int n);

#endif /* __MOLECULE_ASYMROT_EIGSYS_H__ */
