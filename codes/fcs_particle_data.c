/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * fcs_particle_data.c
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#include <config.h>
#include <assert.h>

#include <fcs.h>

#include "fcs_particle_data.h"
#include "data.h"
#include "mathutil.h"
#include "util.h"

/** This function stores all local particles into an array of three consecutive
 *  for coordinates and one array of consecutive doubles for charge and per particle.
 *
 * \note Both arrays should be allocated to the required lengths:
 *  -# \a x: 3 times \a n_local_particles
 *  -# \a q: 1 time \a n_local_particles
 *
 * \note We use fcs_particle_data::n_local_particles for the number to insert.
 *
 * @param P Problem structure reference
 * @param fcs_particles structure with particle arrays
 */
static void
InsertParticlesIntoFCS_ParticleData(struct Problem *P, struct fcs_particle_data *fcs_particles)
{
  int i, j, k;
  struct Particle *p;
  int n_local = 0;

//  fprintf(stdout, "Printing position for fcs of each particle:\n");
  for (i = P->LCS.LowInnerBorder[0]; i < P->LCS.HighInnerBorder[0]; i++)
    for (j = P->LCS.LowInnerBorder[1]; j < P->LCS.HighInnerBorder[1]; j++)
      for (k = P->LCS.LowInnerBorder[2]; k < P->LCS.HighInnerBorder[2]; k++)
        for (p = P->LCS.C[i][j][k]; p; p = p->next) {
          assert (n_local < fcs_particles->n_local_particles);
          // unscaling of positions checked, are again in true [0,P->Box.HMat]^3
          RVec3Mat33(&fcs_particles->positions[n_local * NDIM], p->x, P->Box.HMat);
//          if (MASTER)
//            printf("particle %d: %f %f %f\n", p->Id, fcs_particles->positions[n_local * NDIM+0], fcs_particles->positions[n_local * NDIM+1], fcs_particles->positions[n_local * NDIM+2]);
          fcs_particles->charges[n_local] = p->charge;
          n_local++;
        }
}

/** Allocates memory for arrays contained in fcs_particle_data structure.
 *
 * We use fcs_particles::total_particles as basic number of entries for each array
 *
 * @param fcs_particles pointer to structure
 */
static void
AllocateFCS_ParticleData(struct fcs_particle_data *fcs_particles)
{
  fcs_particles->positions = (fcs_float*)Malloc(3 * fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
  fcs_particles->charges = (fcs_float*)Malloc(fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
  fcs_particles->field = (fcs_float*)Malloc(3 * fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
  fcs_particles->potential = (fcs_float*)Malloc(fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
}

/** Counts the local number of particles.
 *
 * Go through LCStructData and count within cells from LCStructData::LowInnerBorder to
 * LCStructData::HighInnerBorder.
 *
 * @param P Problem structure reference
 * @return number of particles for this process
 */
static int
LocalParticles(struct Problem *P)
{
  int i, j, k;
  struct Particle *p;
  int n_particles = 0;

  for (i = P->LCS.LowInnerBorder[0]; i < P->LCS.HighInnerBorder[0]; i++)
    for (j = P->LCS.LowInnerBorder[1]; j < P->LCS.HighInnerBorder[1]; j++)
      for (k = P->LCS.LowInnerBorder[2]; k < P->LCS.HighInnerBorder[2]; k++)
  for (p = P->LCS.C[i][j][k]; p; p = p->next)
    n_particles++;
  if (verbose)
    fprintf(stderr,"Process %i has %d particles\n", myrank, n_particles);
  return n_particles;
}

/** Initialises the numbers in fcs_particle_data structure for use.
 *
 * @param P pointer to Problem structure
 * @param fcs_particles pointer to structure
 * @param mpi_rank rank of this process
 */
void InitFCS_ParticleData(struct Problem *P, struct fcs_particle_data *fcs_particles, int mpi_rank)
{
  /* count local particles */
  fcs_particles->n_local_particles = LocalParticles(P);

  /* FIXME: this is obviously completely bogus and needs to be fixed */
//  fcs_particles->total_particles = P->Cal.MaxParticleId / P->Par.procs + 10000;
//  fcs_particles->total_particles *= 4;
  MPI_Allreduce( &fcs_particles->n_local_particles, &fcs_particles->total_particles, 1, MPI_INT,
      MPI_SUM, MPI_COMM_WORLD );
  if (MASTER)
    printf( "total_particles = %15d\n", fcs_particles->total_particles );

  MPI_Bcast( &fcs_particles->total_particles, 1, MPI_INT, 0, MPI_COMM_WORLD );

  // get maximum local particle number
  // TODO: this may be set higher, even per process to avoid reallocation in between
  // timesteps when particle number within one process changes
  MPI_Allreduce( &fcs_particles->n_local_particles, &fcs_particles->local_max_particles, 1, MPI_INT,
      MPI_MAX, MPI_COMM_WORLD );

  // allocation
  AllocateFCS_ParticleData(fcs_particles);
}

/** Initialises the numbers in fcs_particle_data structure for use.
 *
 * @param P pointer to Problem structure
 * @param fcs_particles pointer to structure
 * @param mpi_rank rank of this process
 */
void UpdateFCS_ParticleData(struct Problem *P, struct fcs_particle_data *fcs_particles, int mpi_rank)
{
  InsertParticlesIntoFCS_ParticleData(P, fcs_particles);
}

void UpdateFCS_GridData(struct Problem *P, struct fcs_particle_data *fcs_particles)
{
  int n_particles;

  n_particles = LocalParticles(P);

  if(n_particles != fcs_particles->n_local_particles)
    fcs_particles->n_local_particles = n_particles;

  if(n_particles > fcs_particles->local_max_particles){
    fcs_particles->positions = Realloc(fcs_particles->positions, 3 * fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
    fcs_particles->charges = Realloc(fcs_particles->charges, fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
    fcs_particles->field = Realloc(fcs_particles->field, 3 * fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
    fcs_particles->potential = Realloc(fcs_particles->potential, fcs_particles->total_particles * sizeof(fcs_float), "AllocateFCS_ParticleData");
  }
}

/** Frees allocated memory in an fcs_particle_data structure.
 *
 * @param fcs_particles pointer to structure
 */
void FreeFCS_ParticleData(struct fcs_particle_data *fcs_particles)
{
  Free(fcs_particles->positions);
  Free(fcs_particles->charges);
  Free(fcs_particles->field);
  Free(fcs_particles->potential);
}
