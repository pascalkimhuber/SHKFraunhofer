/*
 * FCSData.h
 *
 *  Created on: Jan 31, 2012
 *      Author: heber
 */

#ifndef FCSDATA_H_
#define FCSDATA_H_

#include <fcs.h>

#include "defs.h"
#include "fcs_params.h"

struct fcs_particle_data;
struct Problem;

struct direct_params;
struct ewald_params;
struct fmm_params;
struct memd_params;
struct mmm_params;
struct p2nfft_params;
struct p3m_params;
struct pepc_params;
struct pp3mg_pmg_params;
struct vmg_params;

struct FCSData {
  fcs_int periodic[NDIM];
  fcs_float offset[NDIM];

  //!> cut off radius, this just has to be big enough to encompass any bonddist<=4 contributions
  fcs_float r_cut;
  //!> precision threshold for fcs solvers (most solvers can guarantee this precision for a specific type, see below)
  fcs_float tolerance;
  //!> type of precision for fcs solvers, e.g. absolute or relative, and field, potential, or energy. see FCSDefinitions.h
  fcs_int tolerance_type;

  // internal booleans whether values have been set
  fcs_int r_cut_set;
  fcs_int tolerance_set;
  fcs_int tolerance_type_set;

  //!> FCS solver method to use
  enum FCSMethod method;

  //!> union of all solver-specific parameters
  union {
    struct direct_params *DIRECT;
    struct ewald_params *EWALD;
    struct fmm_params *FMM;
    struct memd_params *MEMD;
    struct mmm_params *MMM;
    struct p2nfft_params *P2NFFT;
    struct p3m_params *P3M;
    struct pepc_params *PEPC;
    struct pp3mg_pmg_params *PP3MG_PMG;
    struct vmg_params *VMG;
  } params;

  //!> particle structure with all internal arrays for the FCS library call
  struct fcs_particle_data *fcs_particles;

  //!> internal handle for FCS library
  FCS fcs_handle;
};

void InitFCSData(struct Problem *P, struct FCSData *FCS);
void FreeFCSData(struct FCSData *FCS);


#endif /* FCSDATA_H_ */
