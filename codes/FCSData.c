/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * FCSData.c
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#include "FCSData.h"

#include "data.h"
#include "defs.h"
#include "fcs_particle_data.h"
#include "util.h"


/** Function to allocate pointers and set initial values for FCSData structure.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData structure to set
 */
void InitFCSData(struct Problem *P, struct FCSData *FCS)
{
  int i;

  /* set periodicity */
  for (i=0; i<NDIM; ++i) {
    FCS->periodic[i] = P->Box.Border[2*i] == periodic;
    FCS->offset[i] = 0;
  }

  FCS->fcs_particles = (struct fcs_particle_data *) Malloc( 1*sizeof(struct FCSData), "InitFCSData" );
}

/** Frees pointers in FCSData structure.
 *
 * @param FCS pointer to FCSData structure to set
 */
void FreeFCSData(struct FCSData *FCS)
{
  FreeFCS_ParticleData(FCS->fcs_particles);
  Free(FCS->fcs_particles);
}
