/*
 * memd_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef MEMD_PARAMS_H_
#define MEMD_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct memd_params
{
	//!> ???
  fcs_int mesh_size;
	//!> ???
  fcs_float lightspeed;
	//!> ???
  fcs_float temperature;
	//!> ???
  fcs_float permittivity;

  // booleans to check whether value was set
  fcs_int mesh_size_set;
  fcs_int lightspeed_set;
  fcs_int temperature_set;
  fcs_int permittivity_set;
};


void InitMEMDParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeMEMDParameters(struct FCSData *FCS);

void parseMEMDParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printMEMDParameters(struct FCSData *FCS, int mpi_rank);
void setMEMDTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* MEMD_PARAMS_H_ */
