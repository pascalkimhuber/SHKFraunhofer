/** @file tremolo-fcs.c
 *  @brief Interface to ScaFaCoS FCS library
 *
 */

#include <config.h>
#include <assert.h>
#include <string.h>
#include <math.h>

/* fcs includes */
#include "fcs.h"
#include "fcs_vmg_p.h"

/* Tremolo includes */
#include "tremolo-fcs.h"
#include "data.h"
#include "errors.h"
#include "parser.h"
#include "coulomb.h"
#include "lcforces.h"
#include "util.h"

#include "FCS/fcs_helpers.h"
#include "FCS/fcs_params.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "FCS/vmg_params.h"

/** Default error handler for FCS module.
 *
 * This function handles errors thrown internally by the FCS module.
 *
 * @param fcs_result result from the module
 */
void HandleFCSError(FCSResult fcs_result)
{
  if (fcs_result) {
    fcsResult_printResult(fcs_result);
    Error(SomeError, "FCS Error");
  }
}

/** Sends the common set of parameters to all other processes from root.
 *
 * @param params pointer to parameter structure
 */
void sendFCSParametersToAll(struct FCSData *params)
{
  MPI_Bcast( &params->periodic[0], 3, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->offset[0], 3, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->r_cut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->tolerance_type, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->r_cut_set, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->tolerance_set, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &params->tolerance_type_set, 1, MPI_INT, 0, MPI_COMM_WORLD );
}


/** Updates Particle::F force array based on FCS result in \a *f.
 *
 * @param P Problem structure reference
 * @param n_local_particles number of particles.
 * @param f forces resulting from FCS module call
 */
static void
UpdateFCSForces(struct Problem *P, int n_local_particles, const fcs_float *f)
{
  int i, j, k;
  struct Particle *p;
  int n_local = 0;

  struct fcs_particle_data *fdata = P->Coul->u.FCS->fcs_particles;

  fprintf(stdout, "Printing Coulomb force of each particle:\n");
  for (i = P->LCS.LowInnerBorder[0]; i < P->LCS.HighInnerBorder[0]; i++)
    for (j = P->LCS.LowInnerBorder[1]; j < P->LCS.HighInnerBorder[1]; j++)
      for (k = P->LCS.LowInnerBorder[2]; k < P->LCS.HighInnerBorder[2]; k++)
	for (p = P->LCS.C[i][j][k]; p; p = p->next) {
	  assert (n_local < n_local_particles);
	  // forces are relative to correct box size, hence don't transform here
	  p->F[0] += f[n_local * NDIM+0]* P->epsilon0inv * p->charge;
	  p->F[1] += f[n_local * NDIM+1]* P->epsilon0inv * p->charge;
	  p->F[2] += f[n_local * NDIM+2]* P->epsilon0inv * p->charge;

	  //update potentials as Scafacos does not use epsilon in the force/potential calculation

	  fdata->potential[n_local] *= 0.5 * P->epsilon0inv * p->charge;

	  n_local++;
	}
}

static double FCSCoulombForce(void *data,
    struct LCForceFunctionData * const LFFD)
{
  struct CoulombData *Coul = data;
  const struct Particle *p = LFFD->p;
  const struct Particle *q = LFFD->q;
  const int bonddist = LFFD->bonddist;
  double F = 0.;
  double G;
  double r_norm = LFFD->r_norm;
  double r_2 = LFFD->r_2;

  // negative of CoulombForce()
  if (bonddist >= 0 && bonddist < 4) {
    G = Coul->P->epsilon0inv * p->charge * q->charge / r_norm;
    F = - G  / r_2;
    if (Coul->calcEnergy)
      Coul->E_k += G;
  }
  return -F;
}


/** Parse FCS entry from the .parameters file in Coulomb section.
 *
 * For now we parse the following parameters:
 * -# cellratio: ration between linked cell and multi grid cells.
 * -# i_degree: interpolation degree of the b-splines
 *
 * For the moment we only support periodic boundary conditions.
 *
 * @param P Problem structure reference
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 * @return 0 - parsing successful.
 */
int ReadFCSCoulombRecord(struct Problem *P, FilePosType *filePos, parse_data *pd)
{
  int i;
  int j;
  struct FCSData *FCS;
  char methodname[256];

  P->Coul->u.FCS = FCS = Malloc(sizeof *FCS, "ReadFCSCoulombRecord");

  // read next entry which gives the method
  FCS->method=0;
  if (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(method):
    		strncpy(methodname, parse_string(filePos, pd), 255);
    	break;
    default:
      parse_error(ERRPARSE, filePos, "ReadFCSCoulombRecord: invalid variable");
      break;
    }
  	for (; FCS->method< fcs_unknown; ++FCS->method) {
  		if (strncmp(methodname, FCSInit[FCS->method].methodname, 255) == 0)
  			break;
  	}
  }
  if (FCS->method != fcs_unknown) {
#ifndef USE_GPL
    if(FCS->method == fcs_p2nfft || FCS->method == fcs_p3m || FCS->method == fcs_mmm || FCS->method == fcs_memd)
      Error(SomeError, "This method is not supported as it is under GPL license.");
#endif
  	FCSInit[FCS->method].parseParamsFn(FCS, filePos, pd);
  } else {
    parse_error(ERRPARSE, filePos, "ReadFCSCoulombRecord: invalid variable");
  }
  for (i=0; i<NDIM; i++) {
    if (P->Box.Border[2*i] == periodic && P->Box.Border[2*i+1] == periodic)
      FCS->periodic[i] = 1;
#if 0
    else if (P->Box.Border[2*i] == reflecting && P->Box.Border[2*i+1] == reflecting)
      FCS->periodic[i] = 0;
#endif
    else
      Error(SomeError, "FCS only supports periodic boundary conditions");
  }
  if (FCS->periodic[0] != FCS->periodic[1] || FCS->periodic[1] != FCS->periodic[2])
    Error(SomeError, "FCS does not work with 1D or 2D periodicity yet");

  P->Cal.CorrectImpulse = 1;

  /* FIXME Array is not initialized if no other short range forces are present, check parameter */

  P->LCS.LFPArray[LCForceParamA] =
      InitLCForceParams(P, P->LCS.LFPArray[LCForceParamA], LCOffset, 14, 1, P->Mol.MaxParticleType, 4);

  for (i=0; i<P->Mol.MaxParticleType; i++)
    for (j=0; j<P->Mol.MaxParticleType; j++)
    {
      /* FIXME define Measure-functions for coloumb and check if everything is thread-safe (last parameter) */
      RegisterLCPotential(P, P->LCS.LFPArray[LCForceParamA], FCSCoulombForce,
	  RemovePotDataDefault, SetMeasureDefault, GetMeasureDefault,
	  i, j, FCS->r_cut, P->Coul, 1, 0);
    }
  return 0;
}

void UpdateFCSStructFactors(struct Problem *P)
{
}

/** Updates the FCS Grid.
 *
 * So far, this is empty.
 *
 * @param P Problem structure reference
 */
void UpdateFCSGrid(struct Problem *P)
{
  UpdateFCS_GridData(P, P->Coul->u.FCS->fcs_particles);
}

/** Initialize FCS parameters.
 *
 * Here, as we now have the particles, we initializes fcs_particle_data and
 * call fcs_common_set which needs the total number of particles.
 *
 * @param P Problem structure reference
 */
void InitFCS(struct Problem *P)
{
  /* result handle */
  FCSResult fcs_result;

  /* MPI variables */
  int mpi_rank;

  struct FCSData *FCS = P->Coul->u.FCS;
  struct fcs_particle_data *fcs_particles = FCS->fcs_particles;

  /* Initializing MPI */
  MPI_Comm_rank( P->Com.comm, &mpi_rank );

  // initialize FCS arrays
  InitFCS_ParticleData(P, FCS->fcs_particles, mpi_rank);

  /* FIXME: UpdateFCSStructFactors and UpdateFCSGrid are not yet implemented but needed for parallel calculations*/
 // if(P->Par.procs > 1 && MASTER)
 //   Error(SomeError, "FCS does not support parallel calculation right now.");

  // set parameters
  printf("fcs_common_set\n");

  fcs_result = fcs_set_common(FCS->fcs_handle, 0,
      P->Box.HMat, &P->Box.HMat[NDIM], &P->Box.HMat[2*NDIM],
      FCS->offset,
      FCS->periodic,
      fcs_particles->total_particles);
  HandleFCSError(fcs_result);
  fcs_set_near_field_flag(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);

  if (FCS->tolerance_set && FCS->tolerance_type_set) {
    printf("Setting tolerance to %.5E of type %d\n", FCS->tolerance, FCS->tolerance_type);
    fcs_result = fcs_set_tolerance(FCS->fcs_handle, FCS->tolerance_type, FCS->tolerance);
  } else if (FCS->tolerance_set) {
    Error(SomeError, "fcs: tolerance set but not what type of tolerance.");
  } else if (FCS->tolerance_type_set) {
    Error(SomeError, "fcs: tolerance type set but not to what value it should be set.");
  }


  /* set and distribute solver specific parameters */
  FCSInit[FCS->method].InitFn(P, FCS, mpi_rank);
}

/** Initiate FCS usage within TREMOLO.
 *
 * \note This comes before InitFCS().
 *
 * We check for cubic box here.
 * We initialize FCSData and the solver-specific parameters. We call fcs_init()
 * and set the dimension via fcs_set_dimension(), finally we call the solver
 * specific setup function.
 *
 * I.e. this is "initialize everything without particles" so far
 *
 * @param P Problem structure reference
 */
void StartFCS(struct Problem *P)
{
  /* result handle */
  FCSResult fcs_result;

  /* MPI variables */
  int mpi_rank;

  struct FCSData *FCS = P->Coul->u.FCS;

  //struct FCSData *FCS = P->Coul->u.FCS;
  if (P->Box.BoxType != BoxIsCube)
    Error(SomeError, "fcs: only implemented for cubic box");

  if (MASTER) {
    fputs("Coulomb potential: using FCS algorithm\n", stderr);
  }

  /* Initializing MPI */
  MPI_Comm_rank( P->Com.comm, &mpi_rank );

  /* set common parameters */
  InitFCSData(P, P->Coul->u.FCS);

  // check FCS method
  assert(FCS->method != fcs_unknown);

  // init with general stuff
  fcs_result = fcs_init(&FCS->fcs_handle, FCSInit[FCS->method].methodname, P->Com.comm);
  HandleFCSError(fcs_result);

  // set dimension
  fcs_result = fcs_set_dimension(FCS->fcs_handle, NDIM);
  HandleFCSError(fcs_result);
}


/** Evaluate forces using the FCS module.
 *
 * Here is the main work horse. We update the particle positions, call fcs_tune and
 * then timed fcs_run() to calculate forces and potential. Subsequently, we use
 * UpdateFCSForces() to set the Particle::F array.
 *
 * @param P Problem structure reference
 */
void FCSForce(struct Problem *P)
{
  /* result handle */
  FCSResult fcs_result;

  /* MPI variables */
  //int mpi_size;
  int mpi_rank;
  double e_sum_local;

  /* Timer variables */
  double t_start = 0, t_stop;

  struct FCSData *FCS = P->Coul->u.FCS;
  struct fcs_particle_data *fcs_particles = FCS->fcs_particles;

  MPI_Comm mpi_comm_cart = P->Com.comm;

  //MPI_Comm_size( P->Com.comm, &mpi_size );
  MPI_Comm_rank( P->Com.comm, &mpi_rank );

  /* ----------------------------------------------------------------------------
   *
   * Updating particles
   * For now, we simply write all particles again into the buffers
   *
   * ---------------------------------------------------------------------------- */

  UpdateFCS_ParticleData(P, FCS->fcs_particles, mpi_rank);

  /* ----------------------------------------------------------------------------
   *
   * Call tune FCS
   *
   * ---------------------------------------------------------------------------- */

  // initiate tuning (inside the module)
  fcs_result = fcs_tune (FCS->fcs_handle,
      fcs_particles->n_local_particles,
      fcs_particles->local_max_particles,
      fcs_particles->positions, fcs_particles->charges);
  HandleFCSError(fcs_result);

  // print active parameters after tuning
  if(verbose)
    FCSInit[FCS->method].printParamsFn(FCS, mpi_rank);

  /* ----------------------------------------------------------------------------
   *
   * Energy and force calculation using FCS
   *
   * ---------------------------------------------------------------------------- */

  /* Starting timer */
  MPI_Barrier( mpi_comm_cart );
  if( mpi_rank == 0 )
    t_start = MPI_Wtime();

  fcs_result = fcs_run(FCS->fcs_handle,
      fcs_particles->n_local_particles, fcs_particles->local_max_particles,
      fcs_particles->positions, fcs_particles->charges, fcs_particles->field, fcs_particles->potential);
  HandleFCSError(fcs_result);

  /* ----------------------------------------------------------------------------
   *
   * Force update
   *
   * ---------------------------------------------------------------------------- */
  UpdateFCSForces(P, fcs_particles->n_local_particles, fcs_particles->field);

  /* Stopping timer */
  MPI_Barrier( mpi_comm_cart );
  if( mpi_rank == 0 ){
    t_stop = MPI_Wtime();
    printf( "Energy and force calculation using %s:              %f s\n", FCSInit[FCS->method].methodname, t_stop-t_start );
  }

  /* ----------------------------------------------------------------------------
   *
   * Energy calculation
   *
   * ---------------------------------------------------------------------------- */
  e_sum_local = sumUpLocalEnergy(fcs_particles);

  // call solver specific setter for total energy
  FCSInit[FCS->method].setTotalEnergyFn(P, FCS, e_sum_local, mpi_rank);
}

/** Cleanup function for the FCS algorithm.
 *
 * We clean up solver-specific parameter structure, FCSData, and also
 * remove the FCS handle via fcs_destroy().
 *
 * @param P pointer to Problem structure
 */
void CleanupFCS(struct Problem *P)
{
  struct FCSData *FCS = P->Coul->u.FCS;

  /* Free memory */
  FCSInit[FCS->method].CleanupFn(FCS);
  FreeFCSData(FCS);

  /* Deallocating storage */
  fcs_destroy(FCS->fcs_handle);
  // Freeing FSCData itself is done in CleanupCoulomb() as ParticleMesh is in same union
}
