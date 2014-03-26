/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * pp3mg_pmg_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "pp3mg_pmg_params.h"

#include "coulomb.h"
#include "data.h"
#include "errors.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "FCS/tremolo-fcs.h"
#include "mathutil.h"
#include "parser.h"
#include "util.h"


/** Sends the set parameters to all other processes from root.
 *
 * @param params pointer to parameter structure
 */
static void
sendPP3MG_PMGParametersToAll(struct pp3mg_pmg_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[8] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[8] = { 1, 1, NDIM, 1, 1, 1, 1, 1 };
  MPI_Aint disp[8];
  MPI_Get_address(&params->degree, &disp[0]);
  MPI_Get_address(&params->ghosts, &disp[1]);
  MPI_Get_address(&params->gridsize[0], &disp[2]);
  MPI_Get_address(&params->max_iterations, &disp[3]);
  MPI_Get_address(&params->degree_set, &disp[4]);
  MPI_Get_address(&params->ghosts_set, &disp[5]);
  MPI_Get_address(&params->gridsize_set, &disp[6]);
  MPI_Get_address(&params->max_iterations_set, &disp[7]);
  MPI_Type_create_struct(8, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates PP3MG_PMG specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitPP3MG_PMGParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct pp3mg_pmg_params *params = FCS->params.PP3MG_PMG;
  FCSResult fcs_result;

  /// send to all other processes
  sendPP3MG_PMGParametersToAll(params);
  sendFCSParametersToAll(FCS);

  if (params->degree_set) {
    fcs_result = fcs_pp3mg_set_degree(FCS->fcs_handle, params->degree);
    HandleFCSError(fcs_result);
  }
  if (params->ghosts_set) {
    fcs_result = fcs_pp3mg_set_ghosts(FCS->fcs_handle, params->ghosts);
    HandleFCSError(fcs_result);
  }
  if (params->gridsize_set) {
    fcs_result = fcs_pp3mg_set_cells_x(FCS->fcs_handle, params->gridsize[0]);
    HandleFCSError(fcs_result);
    fcs_result = fcs_pp3mg_set_cells_y(FCS->fcs_handle, params->gridsize[1]);
    HandleFCSError(fcs_result);
    fcs_result = fcs_pp3mg_set_cells_z(FCS->fcs_handle, params->gridsize[2]);
    HandleFCSError(fcs_result);
  }
  if (params->max_iterations_set) {
    fcs_result = fcs_pp3mg_set_max_iterations(FCS->fcs_handle, params->max_iterations);
    HandleFCSError(fcs_result);
  }
  // TODO: pp3mg does not adhere to fcs_set_tolerance or fcs_set_tolerance_type
  if (FCS->tolerance_set) {
    fcs_result = fcs_pp3mg_set_tol(FCS->fcs_handle, FCS->tolerance);
    HandleFCSError(fcs_result);
  }
  fcs_result = fcs_pp3mg_set_max_particles(FCS->fcs_handle, FCS->fcs_particles->total_particles);
  HandleFCSError(fcs_result);
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}

/** Frees PP3MG_PMG specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreePP3MG_PMGParameters(struct FCSData *FCS)
{
  Free(FCS->params.PP3MG_PMG);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parsePP3MG_PMGParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct pp3mg_pmg_params *params;
  FCS->params.PP3MG_PMG = (struct pp3mg_pmg_params *) Malloc( 1*sizeof(struct pp3mg_pmg_params), "parsePP3MG_PMGParameters");
  params = FCS->params.PP3MG_PMG;

  while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(degree):
      params->degree = parse_int(filePos, pd);
      params->degree_set = 1;
      break;
    case LC_IDENT(ghosts):
      params->ghosts = parse_int(filePos, pd);
      params->ghosts_set = 1;
      break;
    case LC_IDENT(gridsize_x):
      params->gridsize[0] = parse_int(filePos, pd);
      break;
    case LC_IDENT(gridsize_y):
      params->gridsize[1] = parse_int(filePos, pd);
      break;
    case LC_IDENT(gridsize_z):
      params->gridsize[2] = parse_int(filePos, pd);
      break;
    case LC_IDENT(max_iterations):
      params->max_iterations = parse_int(filePos, pd);
      params->max_iterations_set = 1;
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
      break;
    default:
      parse_error(ERRPARSE, filePos, "ReadFCSCoulombRecord: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printPP3MG_PMGParameters(struct FCSData *FCS, int mpi_rank)
{
	struct pp3mg_pmg_params *params = FCS->params.PP3MG_PMG;
  FCSResult fcs_result;
  fcs_float precision;

  // get all params
  fcs_result = fcs_pp3mg_get_degree(FCS->fcs_handle, &params->degree);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_ghosts(FCS->fcs_handle, &params->ghosts);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_cells_x(FCS->fcs_handle, &params->gridsize[0]);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_cells_y(FCS->fcs_handle, &params->gridsize[1]);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_cells_z(FCS->fcs_handle, &params->gridsize[2]);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_max_iterations(FCS->fcs_handle, &params->max_iterations);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pp3mg_get_tol(FCS->fcs_handle, &precision);
  HandleFCSError(fcs_result);
  // TODO: getter is not implemented
//  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
//  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "degree              = %d\n", params->degree );
    printf( "ghosts              = %d\n", params->ghosts );
    printf( "gridsize            = %d %d %d\n", params->gridsize[0], params->gridsize[1], params->gridsize[2] );
    printf( "max_iterations      = %15d\n", params->max_iterations );
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
void setPP3MG_PMGTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct pp3mg_pmg_params *params = FCS->params.PP3MG_PMG;
  int i;
  fcs_float v[9];

  P->Coul->E_l = e_sum_local;

  fcs_get_virial(FCS->fcs_handle, v);

  v[i] *= P->epsilon0inv;


  printf("  virial: %e %e %e %e %e %e %e %e %e\n",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
  // inexpensive error estimator is tr(virial) (see spme.c)
  P->Coul->E_err = v[0] + v[4] + v[8];
}

