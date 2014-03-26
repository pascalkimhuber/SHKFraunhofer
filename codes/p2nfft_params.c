/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * p2nfft_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "p2nfft_params.h"

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
sendP2NFFTParametersToAll(struct p2nfft_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[14] = { MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[14] = { 1, 1, NDIM, 1, 1, NDIM, 1, 1, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[14];
  MPI_Get_address(&params->alpha, &disp[0]);
  MPI_Get_address(&params->epsI, &disp[1]);
  MPI_Get_address(&params->gridsize[0], &disp[2]);
  MPI_Get_address(&params->interpolation_order, &disp[3]);
  MPI_Get_address(&params->m, &disp[4]);
  MPI_Get_address(&params->oversampled_gridsize[0], &disp[5]);
  MPI_Get_address(&params->p, &disp[6]);
  MPI_Get_address(&params->alpha_set, &disp[7]);
  MPI_Get_address(&params->epsI_set, &disp[8]);
  MPI_Get_address(&params->gridsize_set, &disp[9]);
  MPI_Get_address(&params->interpolation_order_set, &disp[10]);
  MPI_Get_address(&params->m_set, &disp[11]);
  MPI_Get_address(&params->oversampled_gridsize_set, &disp[12]);
  MPI_Get_address(&params->p_set, &disp[13]);
  MPI_Type_create_struct(14, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates P2NFFT specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitP2NFFTParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct p2nfft_params *params = FCS->params.P2NFFT;
  FCSResult fcs_result;

  /// send to all other processes
  sendP2NFFTParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (params->alpha_set) {
    fcs_result = fcs_p2nfft_set_alpha(FCS->fcs_handle, params->alpha);
    HandleFCSError(fcs_result);
  }
  if (params->epsI_set) {
    fcs_result = fcs_p2nfft_set_epsI(FCS->fcs_handle, params->epsI);
    HandleFCSError(fcs_result);
  }
  if (params->gridsize_set) {
    fcs_result = fcs_p2nfft_set_grid(FCS->fcs_handle,
        params->gridsize[0], params->gridsize[1], params->gridsize[2]);
    HandleFCSError(fcs_result);
  }
  if (params->interpolation_order_set) {
    fcs_result = fcs_p2nfft_set_pnfft_interpolation_order(FCS->fcs_handle, params->interpolation_order);
    HandleFCSError(fcs_result);
  }
  if (params->m_set) {
    fcs_result = fcs_p2nfft_set_pnfft_m(FCS->fcs_handle, params->m);
    HandleFCSError(fcs_result);
  }
  if (params->oversampled_gridsize_set) {
    fcs_result = fcs_p2nfft_set_oversampled_grid(FCS->fcs_handle,
        params->oversampled_gridsize[0], params->oversampled_gridsize[1], params->oversampled_gridsize[2]);
    HandleFCSError(fcs_result);
  }
  if (params->p_set) {
    fcs_result = fcs_p2nfft_set_p(FCS->fcs_handle, params->p);
    HandleFCSError(fcs_result);
  }
  if (FCS->r_cut_set) {
    fcs_result = fcs_p2nfft_set_r_cut(FCS->fcs_handle, FCS->r_cut);
    HandleFCSError(fcs_result);
  }
  // TODO: P2NFFT does not properly use set_tolerance
  if (FCS->tolerance_set) {
    fcs_result = fcs_p2nfft_set_tolerance(FCS->fcs_handle, FCS->tolerance_type, FCS->tolerance);
    HandleFCSError(fcs_result);
  }
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}

/** Frees P2NFFT specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeP2NFFTParameters(struct FCSData *FCS)
{
  Free(FCS->params.P2NFFT);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseP2NFFTParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct p2nfft_params *params;
  FCS->params.P2NFFT = (struct p2nfft_params *) Malloc( 1*sizeof(struct p2nfft_params), "parseP2NFFTParameters");
  params = FCS->params.P2NFFT;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(alpha):
      params->alpha = parse_double(filePos, pd);
      params->alpha_set = 1;
      break;
    case LC_IDENT(epsI):
      params->epsI = parse_double(filePos, pd);
      params->epsI_set = 1;
      break;
    case LC_IDENT(gridsize_x):
      params->gridsize[0] = parse_int(filePos, pd);
      params->gridsize_set = 1;
      break;
    case LC_IDENT(gridsize_y):
      params->gridsize[1] = parse_int(filePos, pd);
      params->gridsize_set = 1;
      break;
    case LC_IDENT(gridsize_z):
      params->gridsize[2] = parse_int(filePos, pd);
      params->gridsize_set = 1;
      break;
    case LC_IDENT(interpolation_order):
      params->interpolation_order = parse_int(filePos, pd);
      params->interpolation_order_set = 1;
      break;
    case LC_IDENT(m):
      params->m = parse_int(filePos, pd);
      params->m_set = 1;
      break;
    case LC_IDENT(oversampled_gridsize_x):
      params->oversampled_gridsize[0] = parse_int(filePos, pd);
      params->oversampled_gridsize_set = 1;
      break;
    case LC_IDENT(oversampled_gridsize_y):
      params->oversampled_gridsize[1] = parse_int(filePos, pd);
      params->oversampled_gridsize_set = 1;
      break;
    case LC_IDENT(oversampled_gridsize_z):
      params->oversampled_gridsize[2] = parse_int(filePos, pd);
      params->oversampled_gridsize_set = 1;
      break;
    case LC_IDENT(p):
      params->p = parse_int(filePos, pd);
      params->p_set = 1;
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
        parse_error(ERRPARSE, filePos, "parseP2NFFTParameters: tolerance_type not supported.");
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseP2NFFTParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printP2NFFTParameters(struct FCSData *FCS, int mpi_rank)
{
	struct p2nfft_params *params = FCS->params.P2NFFT;
  FCSResult fcs_result;

  // get all params
  fcs_result = fcs_p2nfft_get_alpha(FCS->fcs_handle, &params->alpha);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_epsI(FCS->fcs_handle, &params->epsI);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_grid(FCS->fcs_handle, &params->gridsize[0], &params->gridsize[1], &params->gridsize[2]);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_pnfft_interpolation_order(FCS->fcs_handle, &params->interpolation_order);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_pnfft_m(FCS->fcs_handle, &params->m);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_oversampled_grid(FCS->fcs_handle, &params->oversampled_gridsize[0], &params->oversampled_gridsize[1], &params->oversampled_gridsize[2]);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_p(FCS->fcs_handle, &params->p);
  HandleFCSError(fcs_result);
  fcs_result = fcs_p2nfft_get_r_cut(FCS->fcs_handle, &FCS->r_cut);
  HandleFCSError(fcs_result);
  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "alpha               = %1.5E\n", params->alpha );
    printf( "epsI                = %1.5E\n", params->epsI );
    printf( "gridsize            = %d %d %d\n", params->gridsize[0], params->gridsize[1], params->gridsize[2] );
    printf( "interpolation_order = %d\n", params->interpolation_order );
    printf( "m                   = %d\n", params->m );
    printf( "oversampled_gridsize = %d %d %d\n", params->oversampled_gridsize[0], params->oversampled_gridsize[1], params->oversampled_gridsize[2] );
    printf( "p                   = %d\n", params->p );
    printf( "r_cut               = %.5E\n", FCS->r_cut );
    printf( "tolerance           = %.5E\n", FCS->tolerance );
    // TODO: Tolerance_type is not correctly get with p2nfft?
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
void setP2NFFTTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct p2nfft_params *params = FCS->params.P2NFFT;
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

