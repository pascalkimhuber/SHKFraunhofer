/*
 * ewald_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef EWALD_PARAMS_H_
#define EWALD_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct ewald_params
{
  //!> alpha
  fcs_float alpha;
  //!> kmax
  fcs_int kmax;
  //!> max kmax
  fcs_int maxkmax;

  // booleans to check whether value was set
  fcs_int alpha_set;
  fcs_int kmax_set;
  fcs_int maxkmax_set;
};

void InitEWALDParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeEWALDParameters(struct FCSData *FCS);

void parseEWALDParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printEWALDParameters(struct FCSData *FCS, int mpi_rank);
void setEWALDTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* EWALD_PARAMS_H_ */
