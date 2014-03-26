/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * mmm_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "mmm_params.h"

#include "data.h"
#include "errors.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "parser.h"
#include "util.h"

/** Allocates MMM specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitMMMParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// allocate
  FCS->params.MMM = (struct mmm_params *) Malloc( 1*sizeof(struct mmm_params), "InitMMMParameters");

  /// set to default values
  //struct mmm_params *params = FCS->params.MMM;
  Error(SomeError, "not implemented.");
}

/** Frees MMM specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeMMMParameters(struct FCSData *FCS)
{
  Free(FCS->params.MMM);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseMMMParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  FCS->params.MMM = (struct mmm_params *) Malloc( 1*sizeof(struct mmm_params), "parseMMMParameters");
//	struct mmm_params *params = FCS->params.MMM;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
    case LC_IDENT(degree):
    		break;
    case LC_IDENT(tolerance):
      FCS->tolerance = parse_double(filePos, pd);
      FCS->tolerance_set = 1;
      break;
    case LC_IDENT(tolerance_type):
      FCS->tolerance_type = parse_int(filePos, pd);
      FCS->tolerance_type_set = 1;
      break;
    default:
      parse_error(ERRPARSE, filePos, "parseMMMParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printMMMParameters(struct FCSData *FCS, int mpi_rank)
{
//	struct mmm_params *params = FCS->params.MMM;
  if( mpi_rank == 0 ){
  }
  Error(SomeError, "not implemented.");
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setMMMTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct mmm_params *params = FCS->params.MMM;
  Error(SomeError, "not implemented.");
}

