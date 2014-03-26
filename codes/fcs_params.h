/*
 * fcs_params.h
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#ifndef FCS_PARAMS_H_
#define FCS_PARAMS_H_

#include "parser.h"

enum FCSMethod {
  fcs_direct,
  fcs_ewald,
  fcs_fmm,
  fcs_memd,
  fcs_mmm,
  fcs_pepc,
  fcs_p2nfft,
  fcs_p3m,
  fcs_pp3mg_pmg,
  fcs_vmg,
  fcs_unknown
};

struct FCSData;
struct Problem;

typedef void (InitFCSSolverParamsFn) (struct Problem *P, struct FCSData *FCS, int mpi_rank);
typedef void (parseFCSSolverParamsFn) (struct FCSData *FCS, FilePosType *filePos, parse_data *pd);
typedef void (printFCSSolverParamsFn) (struct FCSData *FCS, int mpi_rank);
typedef void (setTotalEnergyFCSSolverParamsFn) (struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank);
typedef void (FreeFCSSolverParamsFn) (struct FCSData *FCS);

/* order must match enum FCSMethod ! */
typedef struct
{
  enum keywords type;
  /*@observer@*/ const char *name;
  /*@observer@*/ const char *methodname;
  /*@null@*/ InitFCSSolverParamsFn *InitFn;
  /*@null@*/ parseFCSSolverParamsFn *parseParamsFn;
  /*@null@*/ printFCSSolverParamsFn *printParamsFn;
  /*@null@*/ setTotalEnergyFCSSolverParamsFn *setTotalEnergyFn;
  /*@null@*/ FreeFCSSolverParamsFn *CleanupFn;
} FCSInit_array;

extern const FCSInit_array FCSInit[];

#endif /* FCS_PARAMS_H_ */
