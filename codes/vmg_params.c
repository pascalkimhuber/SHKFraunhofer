/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * vmg_params.c
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "coulomb.h"
#include "data.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "FCS/tremolo-fcs.h"
#include "mathutil.h"
#include "parser.h"
#include "util.h"
#include "vmg_params.h"


/** Sends the set parameters to all other processes from root.
 *
 * @param params pointer to parameter structure
 */
static void
sendVMGParametersToAll(struct vmg_params *params)
{
  MPI_Datatype Particletype;
  MPI_Datatype type[16] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklen[16] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[16];
  MPI_Get_address(&params->cycle_type, &disp[0]);
  MPI_Get_address(&params->discretization_order, &disp[1]);
  MPI_Get_address(&params->interpolation_order, &disp[2]);
  MPI_Get_address(&params->max_iterations, &disp[3]);
  MPI_Get_address(&params->max_level, &disp[4]);
  MPI_Get_address(&params->near_field_cells, &disp[5]);
  MPI_Get_address(&params->precision, &disp[6]);
  MPI_Get_address(&params->smoothing_steps, &disp[7]);
  MPI_Get_address(&params->cycle_type_set, &disp[8]);
  MPI_Get_address(&params->discretization_order_set, &disp[9]);
  MPI_Get_address(&params->interpolation_order_set, &disp[10]);
  MPI_Get_address(&params->max_iterations_set, &disp[11]);
  MPI_Get_address(&params->max_level_set, &disp[12]);
  MPI_Get_address(&params->near_field_cells_set, &disp[13]);
  MPI_Get_address(&params->precision_set, &disp[14]);
  MPI_Get_address(&params->smoothing_steps_set, &disp[15]);
  MPI_Type_create_struct(16, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  MPI_Bcast(MPI_BOTTOM, 1, Particletype, 0, MPI_COMM_WORLD );

  MPI_Type_free(&Particletype);
}

/** Allocates VMG specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitVMGParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  struct vmg_params *params = FCS->params.VMG;
  FCSResult fcs_result;

  /// send to all other processes
  sendVMGParametersToAll(params);
  sendFCSParametersToAll(FCS);

  /// initialize FCS method with these parameters
  if (params->cycle_type_set) {
    fcs_result = fcs_vmg_set_cycle_type(FCS->fcs_handle, params->cycle_type);
    HandleFCSError(fcs_result);
  }
  if (params->discretization_order_set) {
    fcs_result = fcs_vmg_set_discretization_order(FCS->fcs_handle, params->discretization_order);
    HandleFCSError(fcs_result);
  }
  if (params->interpolation_order_set) {
    fcs_result = fcs_vmg_set_interpolation_order(FCS->fcs_handle, params->interpolation_order);
    HandleFCSError(fcs_result);
  }
  if (params->max_iterations_set) {
    fcs_result = fcs_vmg_set_max_iterations(FCS->fcs_handle, params->max_iterations);
    HandleFCSError(fcs_result);
  }
  if (params->max_level_set) {
    fcs_result = fcs_vmg_set_max_level(FCS->fcs_handle, params->max_level);
    HandleFCSError(fcs_result);
  }
  if (params->near_field_cells_set) {
    fcs_result = fcs_vmg_set_near_field_cells(FCS->fcs_handle, params->near_field_cells);
    HandleFCSError(fcs_result);
  }
  if (params->precision_set) {
    fcs_result = fcs_vmg_set_precision(FCS->fcs_handle, params->precision);
    HandleFCSError(fcs_result);
  }
  if (params->smoothing_steps_set) {
    fcs_result = fcs_vmg_set_smoothing_steps(FCS->fcs_handle, params->smoothing_steps);
    HandleFCSError(fcs_result);
  }
}

/** Frees VMG specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeVMGParameters(struct FCSData *FCS)
{
  Free(FCS->params.VMG);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseVMGParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  struct vmg_params *params;
  FCS->params.VMG = (struct vmg_params *) Malloc( 1*sizeof(struct vmg_params), "parseVMGParameters");
  params = FCS->params.VMG;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(cycle_type):
      params->cycle_type = parse_int(filePos, pd);
      params->cycle_type_set = 1;
      break;
    case LC_IDENT(discretization_order):
      params->discretization_order = parse_int(filePos, pd);
      params->discretization_order_set = 1;
      break;
    case LC_IDENT(interpolation_order):
      params->interpolation_order = parse_int(filePos, pd);
      params->interpolation_order_set = 1;
      break;
    case LC_IDENT(max_iterations):
      params->max_iterations = parse_int(filePos, pd);
      params->max_iterations_set = 1;
      break;
    case LC_IDENT(max_level):
      params->max_level = parse_int(filePos, pd);
      params->max_level_set = 1;
      break;
    case LC_IDENT(near_field_cells):
      params->near_field_cells = parse_int(filePos, pd);
      params->near_field_cells_set = 1;
      break;
    case LC_IDENT(precision):
      params->precision = parse_double(filePos, pd);
      params->precision_set = 1;
      break;
    case LC_IDENT(r_cut):
      FCS->r_cut = parse_double(filePos, pd);
      FCS->r_cut_set = 1;
      break;
    case LC_IDENT(smoothing_steps):
      params->smoothing_steps = parse_int(filePos, pd);
      params->smoothing_steps_set = 1;
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
      parse_error(ERRPARSE, filePos, "parseVMGParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printVMGParameters(struct FCSData *FCS, int mpi_rank)
{
	struct vmg_params *params = FCS->params.VMG;
  FCSResult fcs_result;

  // get all params
  fcs_result = fcs_vmg_get_cycle_type(FCS->fcs_handle, &params->cycle_type);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_discretization_order(FCS->fcs_handle, &params->discretization_order);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_interpolation_order(FCS->fcs_handle, &params->interpolation_order);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_max_iterations(FCS->fcs_handle, &params->max_iterations);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_max_level(FCS->fcs_handle, &params->max_level);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_near_field_cells(FCS->fcs_handle, &params->near_field_cells);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_precision(FCS->fcs_handle, &params->precision);
  HandleFCSError(fcs_result);
  fcs_result = fcs_vmg_get_smoothing_steps(FCS->fcs_handle, &params->smoothing_steps);
  HandleFCSError(fcs_result);
  // TODO: getter is not implemented
//  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
//  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "cycle_type          = %15d\n", params->cycle_type );
    printf( "discretization_order= %15d\n", params->discretization_order );
    printf( "interpolation degree= %15d\n", params->interpolation_order );
    printf( "max_iterations      = %15d\n", params->max_iterations );
    printf( "max_level           = %15d\n", params->max_level );
    printf( "near_field_cells    = %15d\n", params->near_field_cells );
    printf( "precision           = %.5E\n", params->precision );
    printf( "smoothing_steps     = %15d\n", params->smoothing_steps );
    printf( "tolerance           = %.5E\n", FCS->tolerance );
    printf( "tolerance_type      = %d\n", FCS->tolerance_type );
  }
}

#define minimum3(a,b,c) (MIN(MIN(a,b),c))

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setVMGTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct vmg_params *params = FCS->params.VMG;
  // reduce over all processes
  //MPI_Reduce( (void*) &e_sum_local, (void*) &e_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_cart );

//  /* These values are MPI_Reduce'd later from Tremolo, so give per-cpu ones here */
//  if (mpi_rank == 0)
//    P->Coul->E_k = 1.0/(4.0*M_PI) * 14.0/(5.0*1./(2.0*params->near_field_cells*
//        minimum3(params->grid_size[0],params->grid_size[1],params->grid_size[2])));
//  else
//    P->Coul->E_k = 0.0;
  P->Coul->E_k = 0.;
  P->Coul->E_l = e_sum_local;
  // inexpensive error estimator is tr(virial) (see spme.c)
//  P->Coul->E_err = v[0] + v[4] + v[8];
}
