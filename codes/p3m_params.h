/*
 * p3m_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef P3M_PARAMS_H_
#define P3M_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct p3m_params
{
  //!> ???
  fcs_float alpha;
  //!> ???
  fcs_int cao;
	//!> ???
	fcs_int grid;

  // booleans to check whether value was set
  fcs_int alpha_set;
  fcs_int cao_set;
  fcs_int grid_set;
};

void InitP3MParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeP3MParameters(struct FCSData *FCS);

void parseP3MParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printP3MParameters(struct FCSData *FCS, int mpi_rank);
void setP3MTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* P3M_PARAMS_H_ */
