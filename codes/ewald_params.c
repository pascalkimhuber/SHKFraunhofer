/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * ewald_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "ewald_params.h"

#include "coulomb.h"
#include "data.h"
#include "errors.h"
#include <fcs_ewald_p.h>
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "parser.h"
#include "tremolo-fcs.h"
#include "util.h"

/** Sends the set parameters to all other processes from root.
 *
 * @param params pointer to parameter structure
 */
static void
sendEWALDParametersToAll(struct ewald_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[6] = { MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[6] = { 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[6];
  MPI_Get_address(&params->alpha, &disp[0]);
  MPI_Get_address(&params->kmax, &disp[1]);
  MPI_Get_address(&params->maxkmax, &disp[2]);
  MPI_Get_address(&params->alpha_set, &disp[3]);
  MPI_Get_address(&params->kmax_set, &disp[4]);
  MPI_Get_address(&params->maxkmax, &disp[5]);
  MPI_Type_create_struct(6, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates EWALD specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitEWALDParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct ewald_params *params = FCS->params.EWALD;
  FCSResult fcs_result;

  /// send to all other processes
  sendEWALDParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (params->alpha_set) {
    fcs_result = fcs_ewald_set_alpha(FCS->fcs_handle, params->alpha);
    HandleFCSError(fcs_result);
  }
  if (FCS->r_cut_set) {
    fcs_result = fcs_ewald_set_r_cut(FCS->fcs_handle, FCS->r_cut);
    HandleFCSError(fcs_result);
  }
  if (params->kmax_set) {
    fcs_result = fcs_ewald_set_kmax(FCS->fcs_handle, params->kmax);
    HandleFCSError(fcs_result);
  }
  if (params->maxkmax_set) {
    fcs_result = fcs_ewald_set_maxkmax(FCS->fcs_handle, params->maxkmax);
    HandleFCSError(fcs_result);
  }
}

/** Frees EWALD specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeEWALDParameters(struct FCSData *FCS)
{
  Free(FCS->params.EWALD);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseEWALDParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct ewald_params *params;
  FCS->params.EWALD = (struct ewald_params *) Malloc( 1*sizeof(struct ewald_params), "parseEWALDParameters");
  params = FCS->params.EWALD;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
      case LC_IDENT(alpha):
        params->alpha = parse_double(filePos, pd);
        params->alpha_set = 1;
        break;
    case LC_IDENT(kmax):
      params->kmax = parse_int(filePos, pd);
      params->kmax_set = 1;
      break;
    case LC_IDENT(maxkmax):
      params->maxkmax = parse_int(filePos, pd);
      params->maxkmax_set = 1;
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
      if (FCS->tolerance_type != 5)
        parse_error(ERRPARSE, filePos, "parseEWALDParameters: can only be set to tolerance_type 5");
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseEWALDParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printEWALDParameters(struct FCSData *FCS, int mpi_rank)
{
	struct ewald_params *params = FCS->params.EWALD;
	FCSResult fcs_result;

  // get all params
  fcs_result = fcs_ewald_get_alpha(FCS->fcs_handle, &params->alpha);
  HandleFCSError(fcs_result);
  fcs_result = fcs_ewald_get_r_cut(FCS->fcs_handle, &FCS->r_cut);
  HandleFCSError(fcs_result);
  fcs_result = fcs_ewald_get_kmax(FCS->fcs_handle, &params->kmax);
  HandleFCSError(fcs_result);
  fcs_result = fcs_ewald_get_maxkmax(FCS->fcs_handle, &params->maxkmax);
  HandleFCSError(fcs_result);
  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "alpha               = %f\n", params->alpha);
    printf( "cutoff              = %f\n", FCS->r_cut);
    printf( "kmax                = %d\n", params->kmax);
    printf( "maxkmax             = %d\n", params->maxkmax);
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
void setEWALDTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct ewald_params *params = FCS->params.EWALD;

  P->Coul->E_l = e_sum_local;

  // ewald has not defined any virial functions
//  fcs_float v[9];
//  fcs_ewald_get_virial(FCS->fcs_handle, v);
//  printf("  virial: %e %e %e %e %e %e %e %e %e\n",
//        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
  // inexpensive error estimator is tr(virial) (see spme.c)
//  P->Coul->E_err = v[0] + v[4] + v[8];
}

