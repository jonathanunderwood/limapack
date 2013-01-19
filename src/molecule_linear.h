#ifndef __MOLECULE_LINEAR_H__
#define __MOLECULE_LINEAR_H__

#include "molecule.h"
#include "jmarray.h"
#include "polarizability.h"

typedef struct _linear_molecule
{
  molecule_t parent;
  double B;
  double kT, partfn;
  double coefmin;
  double poptol;
  double oddwt, evenwt;
  polarizability_t * alpha;
  int Jmax;
  int ncoef;
  int nexpval;
  JMarray_int_t * job_status;
} linear_molecule_t;

linear_molecule_t * linear_molecule_ctor();
void linear_molecule_dtor(linear_molecule_t * mol);
void linear_molecule_dispatched_dtor(molecule_t * mol);
molecule_t * linear_molecule_cfg_parse_ctor(const config_setting_t *cfg);

typedef struct _linear_molecule_tdse_worker
{
  int J, M;
} linear_molecule_tdse_worker_t;

#endif /* __MOLECULE_LINEAR_H__ */
