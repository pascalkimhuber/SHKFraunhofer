/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * fcs_helpers.c
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#include <config.h>

#include "data.h"
#include "fcs_helpers.h"
#include "fcs_particle_data.h"

/** Helper method to sum up potential energy over all local particles.
 *
 * @param fcs_particles pointer to fcs_particle_data structure
 * @return energy summed over all particles local to this process
 */
double sumUpLocalEnergy(struct fcs_particle_data *fcs_particles)
{
  double e_sum_local = 0.;
  int p;

  // sum up local energy
  for( p = 0; p < fcs_particles->n_local_particles; p++ )
    e_sum_local += fcs_particles->potential[p];

  return e_sum_local;
}
