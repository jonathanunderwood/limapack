#ifndef __MOLECULE_LINEAR_H__
#define __MOLECULE_LINEAR_H__

#include "molecule.h"

typedef struct _linear_molecule linear_molecule_t;

linear_molecule_t * linear_molecule_ctor();
void linear_molecule_dtor(molecule_t * mol);
molecule_t * linear_molecule_cfg_parse_ctor(const config_setting_t *cfg);

#endif /* __MOLECULE_LINEAR_H__ */
