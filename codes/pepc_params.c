/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * pepc_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "pepc_params.h"

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
sendP2NFFTParametersToAll(struct pepc_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[14] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[14] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[14];
  MPI_Get_address(&params->debuglevel, &disp[0]);
  MPI_Get_address(&params->dipole_correction, &disp[1]);
  MPI_Get_address(&params->epsilon, &disp[2]);
  MPI_Get_address(&params->load_balancing, &disp[3]);
  MPI_Get_address(&params->npm, &disp[4]);
  MPI_Get_address(&params->num_walk_threads, &disp[5]);
  MPI_Get_address(&params->theta, &disp[6]);
  MPI_Get_address(&params->debuglevel_set, &disp[7]);
  MPI_Get_address(&params->dipole_correction_set, &disp[8]);
  MPI_Get_address(&params->epsilon_set, &disp[9]);
  MPI_Get_address(&params->load_balancing_set, &disp[10]);
  MPI_Get_address(&params->npm_set, &disp[11]);
  MPI_Get_address(&params->num_walk_threads_set, &disp[12]);
  MPI_Get_address(&params->theta_set, &disp[13]);
  MPI_Type_create_struct(14, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates PEPC specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitPEPCParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
  struct pepc_params *params = FCS->params.PEPC;
  FCSResult fcs_result;

  /// send to all other processes
  sendP2NFFTParametersToAll(params);
  sendFCSParametersToAll(FCS);

  // set debuglevel is not implemented
//  if (params->debuglevel_set) {
//    fcs_result = fcs_pepc_set_debuglevel(FCS->fcs_handle, params->debuglevel);
//    HandleFCSError(fcs_result);
//  }
  if (params->dipole_correction_set) {
    fcs_result = fcs_pepc_set_dipole_correction(FCS->fcs_handle, params->dipole_correction);
    HandleFCSError(fcs_result);
  }
  if (params->epsilon_set) {
    fcs_result = fcs_pepc_set_epsilon(FCS->fcs_handle, params->epsilon);
    HandleFCSError(fcs_result);
  }
  if (params->load_balancing_set) {
    fcs_result = fcs_pepc_set_load_balancing(FCS->fcs_handle, params->load_balancing);
    HandleFCSError(fcs_result);
  }
  if (params->npm_set) {
    fcs_result = fcs_pepc_set_npm(FCS->fcs_handle, params->npm);
    HandleFCSError(fcs_result);
  }
  if (params->num_walk_threads) {
    fcs_result = fcs_pepc_set_num_walk_threads(FCS->fcs_handle, params->num_walk_threads);
    HandleFCSError(fcs_result);
  }
  if (params->theta_set) {
    fcs_result = fcs_pepc_set_theta(FCS->fcs_handle, params->theta);
    HandleFCSError(fcs_result);
  }
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}
/** Frees PEPC specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreePEPCParameters(struct FCSData *FCS)
{
  Free(FCS->params.PEPC);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parsePEPCParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct pepc_params *params;
  FCS->params.PEPC = (struct pepc_params *) Malloc( 1*sizeof(struct pepc_params), "parsePEPCParameters");
  params = FCS->params.PEPC;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(debuglevel):
      params->debuglevel = parse_int(filePos, pd);
      params->debuglevel_set = 1;
      break;
    case LC_IDENT(dipole_correction):
      params->dipole_correction = parse_int(filePos, pd);
      params->dipole_correction_set = 1;
      break;
    case LC_IDENT(epsilon):
      params->epsilon = parse_double(filePos, pd);
      params->epsilon_set = 1;
      break;
    case LC_IDENT(load_balancing):
      params->load_balancing = parse_int(filePos, pd);
      params->load_balancing_set = 1;
      break;
    case LC_IDENT(npm):
      params->npm = parse_double(filePos, pd);
      params->npm_set = 1;
      break;
    case LC_IDENT(num_walk_threads):
      params->num_walk_threads = parse_double(filePos, pd);
      params->num_walk_threads_set = 1;
      break;
    case LC_IDENT(r_cut):
      FCS->r_cut = parse_double(filePos, pd);
      FCS->r_cut_set = 1;
      break;
    case LC_IDENT(theta):
      params->theta = parse_double(filePos, pd);
      params->theta_set = 1;
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
      parse_error(ERRPARSE, filePos, "parsePEPCParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printPEPCParameters(struct FCSData *FCS, int mpi_rank)
{
	struct pepc_params *params = FCS->params.PEPC;
  FCSResult fcs_result;

  // get all params
//  fcs_result = fcs_pepc_get_debuglevel(FCS->fcs_handle, &params->debuglevel);
//  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_dipole_correction(FCS->fcs_handle, &params->dipole_correction);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_epsilon(FCS->fcs_handle, &params->epsilon);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_load_balancing(FCS->fcs_handle, &params->load_balancing);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_npm(FCS->fcs_handle, &params->npm);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_num_walk_threads(FCS->fcs_handle, &params->num_walk_threads);
  HandleFCSError(fcs_result);
  fcs_result = fcs_pepc_get_theta(FCS->fcs_handle, &params->theta);
  HandleFCSError(fcs_result);
  // TODO: getter is not implemented
//  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
//  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
//    printf( "debuglevel          = %d\n", params->debuglevel );
    printf( "dipole_correction   = %d\n", params->dipole_correction );
    printf( "epsilon             = %.5E\n", params->epsilon );
    printf( "load_balancing      = %d\n", params->load_balancing );
    printf( "npm                 = %.5E\n", params->npm );
    printf( "num_walk_threads    = %d\n", params->num_walk_threads );
    printf( "theta               = %.5E\n", params->theta ) ;
    printf( "tolerance           = unknown\n" );// = .5E\n", FCS->tolerance );
    printf( "tolerance_type      = unknown\n" );// = .5E\n", FCS->tolerance_type );
  }
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setPEPCTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct pepc_params *params = FCS->params.PEPC;
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

