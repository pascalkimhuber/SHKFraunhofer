/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * fmm_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "fmm_params.h"

#include "coulomb.h"
#include "data.h"
#include "errors.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "FCS/tremolo-fcs.h"
#include "parser.h"
#include "util.h"


/** Sends the set parameters to all other processes from root.
 *
 * @param params pointer to parameter structure
 */
static void
sendFMMParametersToAll(struct fmm_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[12] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[12] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[12];
  MPI_Get_address(&params->balanceload, &disp[0]);
  MPI_Get_address(&params->dipole_correction, &disp[1]);
  MPI_Get_address(&params->maxdepth, &disp[2]);
  MPI_Get_address(&params->potential, &disp[3]);
  MPI_Get_address(&params->radius, &disp[4]);
  MPI_Get_address(&params->unroll_limit, &disp[5]);
  MPI_Get_address(&params->balanceload_set, &disp[6]);
  MPI_Get_address(&params->dipole_correction_set, &disp[7]);
  MPI_Get_address(&params->maxdepth_set, &disp[8]);
  MPI_Get_address(&params->potential_set, &disp[9]);
  MPI_Get_address(&params->radius_set, &disp[10]);
  MPI_Get_address(&params->unroll_limit_set, &disp[11]);
  MPI_Type_create_struct(12, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates FMM specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitFMMParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct fmm_params *params = FCS->params.FMM;
  FCSResult fcs_result;

  /// send to all other processes
  sendFMMParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (params->balanceload_set) {
    fcs_result = fcs_fmm_set_balanceload(FCS->fcs_handle, params->balanceload);
    HandleFCSError(fcs_result);
  }
  if (params->dipole_correction_set) {
    fcs_result = fcs_fmm_set_dipole_correction(FCS->fcs_handle, params->dipole_correction);
    HandleFCSError(fcs_result);
  }
  if (params->maxdepth_set) {
    fcs_result = fcs_fmm_set_maxdepth(FCS->fcs_handle, params->maxdepth);
    HandleFCSError(fcs_result);
  }
  if (params->potential_set) {
    fcs_result = fcs_fmm_set_potential(FCS->fcs_handle, params->potential);
    HandleFCSError(fcs_result);
  }
  if (params->radius_set) {
    fcs_result = fcs_fmm_set_cusp_radius(FCS->fcs_handle, params->radius);
    HandleFCSError(fcs_result);
  }
  if (params->unroll_limit_set) {
    fcs_result = fcs_fmm_set_unroll_limit(FCS->fcs_handle, params->unroll_limit);
    HandleFCSError(fcs_result);
  }
  if (FCS->tolerance_type_set) {
    switch (FCS->tolerance_type_set) {
      // TODO: Somehow absrel is always 2 and deltaE is ignored with FMM
      case FCS_TOLERANCE_TYPE_ENERGY:
        fcs_result = fcs_fmm_set_tolerance_energy(FCS->fcs_handle, FCS->tolerance);
        HandleFCSError(fcs_result);
        fcs_result = fcs_fmm_set_absrel(FCS->fcs_handle, FCS->tolerance_type);
        HandleFCSError(fcs_result);
        break;
      case FCS_TOLERANCE_TYPE_ENERGY_REL:
        fcs_result = fcs_fmm_set_absrel(FCS->fcs_handle, FCS->tolerance_type);
        HandleFCSError(fcs_result);
        break;
      default:
        break;
    }
  }
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}

/** Frees FMM specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeFMMParameters(struct FCSData *FCS)
{
  Free(FCS->params.FMM);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseFMMParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct fmm_params *params;
  FCS->params.FMM = (struct fmm_params *) Malloc( 1*sizeof(struct fmm_params), "parseFMMParameters");
  params = FCS->params.FMM;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(balanceload):
      params->balanceload = parse_int(filePos, pd);
      params->balanceload_set = 1;
      break;
    case LC_IDENT(dipole_correction):
      params->dipole_correction = parse_int(filePos, pd);
      params->dipole_correction_set = 1;
      break;
    case LC_IDENT(maxdepth):
      params->maxdepth = parse_int(filePos, pd);
      params->maxdepth_set = 1;
      break;
    case LC_IDENT(potential):
      params->potential = parse_int(filePos, pd);
      params->potential_set = 1;
      break;
    case LC_IDENT(r_cut):
      FCS->r_cut = parse_double(filePos, pd);
      FCS->r_cut_set = 1;
      break;
    case LC_IDENT(radius):
      params->radius = parse_double(filePos, pd);
      params->radius_set = 1;
      break;
    case LC_IDENT(tolerance):
      FCS->tolerance = parse_double(filePos, pd);
      FCS->tolerance_set = 1;
      break;
    case LC_IDENT(tolerance_type):
      FCS->tolerance_type = parse_int(filePos, pd);
      FCS->tolerance_type_set = 1;
      if ((FCS->tolerance_type != FCS_TOLERANCE_TYPE_ENERGY)
          && (FCS->tolerance_type != FCS_TOLERANCE_TYPE_ENERGY_REL))
        parse_error(ERRPARSE, filePos,
            "parseFMMParameters: tolerance type can only be set to absolute (1) or relative (2) energy error");
      break;
    case LC_IDENT(unroll_limit):
      params->unroll_limit = parse_int(filePos, pd);
      params->unroll_limit_set = 1;
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseFMMParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printFMMParameters(struct FCSData *FCS, int mpi_rank)
{
	struct fmm_params *params = FCS->params.FMM;
	FCSResult fcs_result;
	// internal values for tolerance
	fcs_int absrel;
	fcs_float deltaE;

  // get all params
  fcs_result = fcs_fmm_get_absrel(FCS->fcs_handle, &absrel);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_balanceload(FCS->fcs_handle, &params->balanceload);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_cusp_radius(FCS->fcs_handle, &params->radius);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_tolerance_energy(FCS->fcs_handle, &deltaE);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_dipole_correction(FCS->fcs_handle, &params->dipole_correction);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_maxdepth(FCS->fcs_handle, &params->maxdepth);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_potential(FCS->fcs_handle, &params->potential);
  HandleFCSError(fcs_result);
  fcs_result = fcs_fmm_get_unroll_limit(FCS->fcs_handle, &params->unroll_limit);
  HandleFCSError(fcs_result);
  // TODO: getter is not implemented
//  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
//  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "absrel              = %15d\n", absrel );
    printf( "balanceload         = %lld\n", params->balanceload );
    printf( "deltaE              = %.5E\n", deltaE );
    printf( "dipole_correction   = %d\n", params->dipole_correction );
    printf( "maxdepth            = %lld\n", params->maxdepth );
    printf( "potential type      = %d\n", params->potential );
    printf( "radius              = %.5E\n", params->radius );
    if (FCS->tolerance_set)
      printf( "tolerance           = %.5E\n", FCS->tolerance );
    if (FCS->tolerance_type_set)
      printf( "tolerance_type      = %d\n", FCS->tolerance_type );
    printf( "unroll_limit        = %lld\n", params->unroll_limit );
  }
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setFMMTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct fmm_params *params = FCS->params.FMM;
  int i;
  fcs_float v[9];

  P->Coul->E_l = e_sum_local;

  fcs_get_virial(FCS->fcs_handle, v);

  for(i = 0; i < 9; i++)
    v[i] *= P->epsilon0inv;

  printf("  virial: %e %e %e %e %e %e %e %e %e\n",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
  // inexpensive error estimator is tr(virial) (see spme.c)
  P->Coul->E_err = v[0] + v[4] + v[8];
}

