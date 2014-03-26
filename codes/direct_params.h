/*
 * direct_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef DIRECT_PARAMS_H_
#define DIRECT_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct direct_params
{
	//!> number of periodic images used in each (periodic) dimension
	fcs_int periodic_images[NDIM];

	// booleans to check whether value was set
  fcs_int periodic_images_set;
};

void InitDIRECTParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeDIRECTParameters(struct FCSData *FCS);

void parseDIRECTParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printDIRECTParameters(struct FCSData *FCS, int mpi_rank);
void setDIRECTTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* DIRECT_PARAMS_H_ */
