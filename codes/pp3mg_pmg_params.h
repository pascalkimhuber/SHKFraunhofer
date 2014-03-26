/*
 * pp3mg_pmg_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef PP3MG_PMG_PARAMS_H_
#define PP3MG_PMG_PARAMS_H_

#include "data.h"
#include "parser.h"

struct FCSData;
struct Problem;

struct pp3mg_pmg_params
{
  //!> interpolation degree
  fcs_int degree;
  //!> number of ghost cells
  fcs_int ghosts;
  //!> number of grid cells per axis
  fcs_int gridsize[NDIM];
  //!> maximum number of iterations
  fcs_int max_iterations;

  // booleans to check whether value was set
  fcs_int degree_set;
  fcs_int ghosts_set;
  fcs_int gridsize_set;
  fcs_int max_iterations_set;
};

void InitPP3MG_PMGParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreePP3MG_PMGParameters(struct FCSData *FCS);

void parsePP3MG_PMGParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printPP3MG_PMGParameters(struct FCSData *FCS, int mpi_rank);
void setPP3MG_PMGTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* PP3MG_PMG_PARAMS_H_ */
