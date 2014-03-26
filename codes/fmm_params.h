/*
 * fmm_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef FMM_PARAMS_H_
#define FMM_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct fmm_params
{
  //!> load balancing routine (0/1)
  long long balanceload;
	//!> either set to FCS_FMM_NO_DIPOLE_CORRECTION, FCS_FMM_STANDARD_DIPOLE_CORRECTION, or FCS_FMM_ACTIVE_DIPOLE_CORRECTION
	fcs_int dipole_correction;
  //!> maximum tree depth for FMM (0-19)
  long long maxdepth;
	//!> either set to FCS_FMM_COULOMB or FCS_FMM_CUSP
	fcs_int potential;
	//!> radius for cusp potential
	fcs_float radius;
  //!> limit for unrolled functions (0-50)
  long long unroll_limit;

  // booleans to check whether value was set
  fcs_int balanceload_set;
  fcs_int dipole_correction_set;
  fcs_int maxdepth_set;
  fcs_int potential_set;
  fcs_int radius_set;
  fcs_int unroll_limit_set;
};

void InitFMMParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeFMMParameters(struct FCSData *FCS);

void parseFMMParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printFMMParameters(struct FCSData *FCS, int mpi_rank);
void setFMMTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* FMM_PARAMS_H_ */
