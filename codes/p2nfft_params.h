/*
 * p2nfft_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef P2NFFT_PARAMS_H_
#define P2NFFT_PARAMS_H_

#include "parser.h"

struct FCSData;
struct Problem;

struct p2nfft_params
{
  //!> alpha ???
  fcs_float alpha;
  //!> epsI ???
  fcs_float epsI;
  //!> size of fft grid
  fcs_int gridsize[NDIM];
  //!> interpolation order
  fcs_int interpolation_order;
	//!> m ???
	fcs_int m;
  //!> size of oversampled fft grid ???
  fcs_int oversampled_gridsize[NDIM];
	//!> p ???
	fcs_int p;

  // booleans to check whether value was set
  fcs_int alpha_set;
  fcs_int epsI_set;
  fcs_int gridsize_set;
  fcs_int interpolation_order_set;
  fcs_int m_set;
  fcs_int oversampled_gridsize_set;
  fcs_int p_set;
};


void InitP2NFFTParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank);
void FreeP2NFFTParameters(struct FCSData *FCS);

void parseP2NFFTParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
void printP2NFFTParameters(struct FCSData *FCS, int mpi_rank);
void setP2NFFTTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);

#endif /* P2NFFT_PARAMS_H_ */
