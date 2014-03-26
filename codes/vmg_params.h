/*
 * \file vmg_params.h
 *
 * This file contains structures and definitions for the VMG method within the
 * ScaFaCoS library.
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#ifndef VMG_PARAMS_H_
#define VMG_PARAMS_H_

#include "defs.h"
#include "parser.h"

struct FCSData;
struct Problem;

//!> struct with all available parameters of the method
struct vmg_params {
  //!> multi-grid cycle type specifier, e.g. 1 = V-cycle, 2 = W-cycle, ...
  fcs_int cycle_type;
  //!> order of discretization scheme for partial differential equation solver
  fcs_int discretization_order;
  //!> Degree of interpolation polynomial
  fcs_int interpolation_order;
  //!> maximum number of solver iterations
  fcs_int max_iterations;
  //!> maximum number of multi-grid levels
  fcs_int max_level;
  //!> range of smearing of charges over near field cells, must remain with this and at most neighboring process
  fcs_int near_field_cells;
  //!> threshold for (absolute and relative) resdiual for solver iteration
  fcs_float precision;
  //!> number of pre- and post-smoothing steps
  fcs_int smoothing_steps;

  // booleans to check whether value was set
  fcs_int cycle_type_set;
  fcs_int discretization_order_set;
  fcs_int interpolation_order_set;
  fcs_int max_iterations_set;
  fcs_int max_level_set;
  fcs_int near_field_cells_set;
  fcs_int precision_set;
  fcs_int smoothing_steps_set;
};

void InitVMGParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeVMGParameters(struct FCSData *FCS);

void parseVMGParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printVMGParameters(struct FCSData *FCS, int mpi_rank);
void setVMGTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* VMG_PARAMS_H_ */
