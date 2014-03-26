/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * p3m_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "p3m_params.h"

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
sendP2NFFTParametersToAll(struct p3m_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[6] = { MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[6] = { 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[6];
  MPI_Get_address(&params->alpha, &disp[0]);
  MPI_Get_address(&params->cao, &disp[1]);
  MPI_Get_address(&params->grid, &disp[2]);
  MPI_Get_address(&params->alpha_set, &disp[3]);
  MPI_Get_address(&params->cao_set, &disp[4]);
  MPI_Get_address(&params->grid_set, &disp[5]);
  MPI_Type_create_struct(6, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates P3M specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitP3MParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct p3m_params *params = FCS->params.P3M;
  FCSResult fcs_result;

  /// send to all other processes
  sendP2NFFTParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (params->alpha_set) {
    fcs_result = fcs_p3m_set_alpha(FCS->fcs_handle, params->alpha);
    HandleFCSError(fcs_result);
  }
  if (params->cao_set) {
    fcs_p3m_set_cao(FCS->fcs_handle, params->cao);
    HandleFCSError(fcs_result);
  }
  if (params->grid_set) {
    fcs_p3m_set_grid(FCS->fcs_handle, params->grid);
    HandleFCSError(fcs_result);
  }
  if (FCS->r_cut_set) {
    fcs_result = fcs_p3m_set_r_cut(FCS->fcs_handle, FCS->r_cut);
    HandleFCSError(fcs_result);
  }
  // TODO: P3M does not adhere fcs_set_tolerance? crash when not set with internal method
  if (FCS->tolerance_set) {
    fcs_result = fcs_p3m_set_tolerance_field(FCS->fcs_handle, FCS->tolerance);
    HandleFCSError(fcs_result);
  }

  // we need total energy and virial
  fcs_p3m_require_total_energy(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}

/** Frees P3M specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeP3MParameters(struct FCSData *FCS)
{
  Free(FCS->params.P3M);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseP3MParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  FCS->params.P3M = (struct p3m_params *) Malloc( 1*sizeof(struct p3m_params), "parseP3MParameters");
	struct p3m_params *params = FCS->params.P3M;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(alpha):
      params->alpha = parse_double(filePos, pd);
      params->alpha_set = 1;
      break;
    case LC_IDENT(cao):
      params->cao = parse_int(filePos, pd);
      params->cao_set = 1;
      break;
    case LC_IDENT(grid):
      params->grid = parse_int(filePos, pd);
      params->grid_set = 1;
      break;
    case LC_IDENT(r_cut):
      FCS->r_cut = parse_double(filePos, pd);
      FCS->r_cut_set = 1;
      break;
    case LC_IDENT(tolerance):
      FCS->tolerance = parse_double(filePos, pd);
      FCS->tolerance_set = 1;
      break;
    case LC_IDENT(tolerance_type):
      FCS->tolerance_type = parse_int(filePos, pd);
      FCS->tolerance_type_set = 1;
      if (FCS->tolerance_type != FCS_TOLERANCE_TYPE_FIELD)
        parse_error(ERRPARSE, filePos, "parseP3MParameters: only absolute field tolerance_type supported.");
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseP3MParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printP3MParameters(struct FCSData *FCS, int mpi_rank)
{
	struct p3m_params *params = FCS->params.P3M;
  FCSResult fcs_result;
  fcs_float precision;

  // get all params
  fcs_result = fcs_p3m_get_alpha(FCS->fcs_handle, &params->alpha);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p3m_get_cao(FCS->fcs_handle, &params->cao);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p3m_get_grid(FCS->fcs_handle, &params->grid);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p3m_get_tolerance_field(FCS->fcs_handle, &precision);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p3m_get_r_cut(FCS->fcs_handle, &FCS->r_cut);
  HandleFCSError(fcs_result);
  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "alpha               = %.5E\n", params->alpha );
    printf( "cao                 = %15d\n", params->cao );
    printf( "grid                = %15d\n", params->grid );
    printf( "precision           = %.5E\n", precision );
    printf( "r_cut               = %.5E\n", FCS->r_cut );
    printf( "tolerance           = %.5E\n", FCS->tolerance );
    printf( "tolerance_type      = %d\n", FCS->tolerance_type );
  }
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setP3MTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct p3m_params *params = FCS->params.P3M;
  int i;

  P->Coul->E_l = e_sum_local;

  fcs_float v[9];
  fcs_get_virial(FCS->fcs_handle, v);

  v[i] *= P->epsilon0inv;


  printf("  virial: %e %e %e %e %e %e %e %e %e\n",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
  // inexpensive error estimator is tr(virial) (see spme.c)
  P->Coul->E_err = v[0] + v[4] + v[8];
}


