/*
 * pepc_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef PEPC_PARAMS_H_
#define PEPC_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct pepc_params
{
  //!> debuglevel
  fcs_int debuglevel;
  //!> dipole_correction switch ???
  fcs_int dipole_correction;
	//!> epsilon ???
	fcs_float epsilon;
	//!> load balancing switch, only when no particle reordering is done elsewhere
	fcs_int load_balancing;
	//!> npn ???
	fcs_float npm;
	//!> number of walk threads per MPI rank
	fcs_int num_walk_threads;
	//!> ???
	fcs_float theta;

  // booleans to check whether value was set
  fcs_int debuglevel_set;
  fcs_int dipole_correction_set;
  fcs_int epsilon_set;
  fcs_int load_balancing_set;
  fcs_int npm_set;
  fcs_int num_walk_threads_set;
  fcs_int theta_set;
};

void InitPEPCParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreePEPCParameters(struct FCSData *FCS);

void parsePEPCParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printPEPCParameters(struct FCSData *FCS, int mpi_rank);
void setPEPCTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);


#endif /* PEPC_PARAMS_H_ */
