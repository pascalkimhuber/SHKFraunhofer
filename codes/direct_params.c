/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * direct_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "direct_params.h"

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
sendDIRECTParametersToAll(struct direct_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[2] = { MPI_INT, MPI_INT };
  int blocklen[2] = { NDIM, 1 };
  MPI_Aint disp[2];
  MPI_Get_address(&params->periodic_images, &disp[0]);
  MPI_Get_address(&params->periodic_images_set, &disp[1]);
  MPI_Type_create_struct(2, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates DIRECT specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitDIRECTParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct direct_params *params = FCS->params.DIRECT;
  FCSResult fcs_result;

  /// send to all other processes
  sendDIRECTParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (FCS->r_cut_set) {
    fcs_result = fcs_direct_set_cutoff(FCS->fcs_handle, FCS->r_cut);
    HandleFCSError(fcs_result);
  }
  // TODO: set_periodic_images not implemented
  if (params->periodic_images_set) {
    fcs_result = fcs_direct_set_periodic_images(FCS->fcs_handle, &params->periodic_images[0]);
    HandleFCSError(fcs_result);
  }
}

/** Frees DIRECT specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeDIRECTParameters(struct FCSData *FCS)
{
  Free(FCS->params.DIRECT);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseDIRECTParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct direct_params *params;

  FCS->params.DIRECT = (struct direct_params *) Malloc( 1*sizeof(struct direct_params), "parseDIRECTParameters");
  params = FCS->params.DIRECT;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(periodic_images_x):
      params->periodic_images[0] = parse_int(filePos, pd);
      params->periodic_images_set = 1;
      break;
    case LC_IDENT(periodic_images_y):
      params->periodic_images[1] = parse_int(filePos, pd);
      params->periodic_images_set = 1;
      break;
    case LC_IDENT(periodic_images_z):
      params->periodic_images[2] = parse_int(filePos, pd);
      params->periodic_images_set = 1;
      break;
    case LC_IDENT(r_cut):
      FCS->r_cut = parse_double(filePos, pd);
      FCS->r_cut_set = 1;
      break;
    case LC_IDENT(tolerance):
      parse_error(ERRPARSE, filePos, "parseDIRECTParameters: direct solver has always machine precision");
      break;
    case LC_IDENT(tolerance_type):
      parse_error(ERRPARSE, filePos, "parseDIRECTParameters: direct solver has always machine precision");
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseDIRECTParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printDIRECTParameters(struct FCSData *FCS, int mpi_rank)
{
	struct direct_params *params = FCS->params.DIRECT;
	FCSResult fcs_result;

	// get all params
	fcs_result = fcs_direct_get_cutoff(FCS->fcs_handle, &FCS->r_cut);
  HandleFCSError(fcs_result);
  fcs_result = fcs_direct_get_periodic_images(FCS->fcs_handle, &params->periodic_images[0]);
  HandleFCSError(fcs_result);
  // TODO: getter is not implemented
  //  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
  //  HandleFCSError(fcs_result);

  // print params
  if (mpi_rank == 0)
  {
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "cutoff              = %f\n", FCS->r_cut);
    printf( "periodic_images     = %d %d %d\n", params->periodic_images[0], params->periodic_images[1], params->periodic_images[2]);
    printf( "tolerance           = machine precision\n" );
  }
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setDIRECTTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct direct_params *params = FCS->params.DIRECT;

  P->Coul->E_l = e_sum_local;
}

