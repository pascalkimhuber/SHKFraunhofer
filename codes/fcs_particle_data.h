/*
 * fcs_particle_data.h
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#ifndef FCS_PARTICLE_DATA_H_
#define FCS_PARTICLE_DATA_H_

#include <fcs.h>

struct Problem;
struct fcs_particle_data;

//!> Particle data
struct fcs_particle_data {
  //!> linear array of 3 values of coordinates of each atom
  fcs_float *positions;
  //!> linear array of charge per atom
  fcs_float *charges;
  //!> linear array of 3 values of force vector per atom
  fcs_float *field;
  //!> linear array of potential at each atom
  fcs_float *potential;
  //!> local number of particles
  int n_local_particles;
  //!> local maximum(!) number of particles
  int local_max_particles;
  //!> total number of all particles
  int total_particles;
};

void InitFCS_ParticleData(struct Problem *P, struct fcs_particle_data *fcs_particles, int mpi_rank);
void UpdateFCS_ParticleData(struct Problem *P, struct fcs_particle_data *fcs_particles, int mpi_rank);
void UpdateFCS_GridData(struct Problem *P, struct fcs_particle_data *fcs_particles);
void FreeFCS_ParticleData(struct fcs_particle_data *fcs_particles);


#endif /* FCS_PARTICLE_DATA_H_ */
