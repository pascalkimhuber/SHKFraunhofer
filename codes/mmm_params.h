/*
 * mmm_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef MMM_PARAMS_H_
#define MMM_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct mmm_params
{

};

void InitMMMParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeMMMParameters(struct FCSData *FCS);

void parseMMMParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printMMMParameters(struct FCSData *FCS, int mpi_rank);
void setMMMTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* MMM_PARAMS_H_ */
