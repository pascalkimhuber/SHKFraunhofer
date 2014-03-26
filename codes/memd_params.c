/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * memd_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include <config.h>

#include <fcs.h>

#include "memd_params.h"

#include "coulomb.h"
#include "data.h"
#include "errors.h"
#include "FCS/fcs_particle_data.h"
#include "FCS/FCSData.h"
#include "FCS/tremolo-fcs.h"
#include "parser.h"
#include "util.h"

/** Allocates MEMD specific memory inside FCSData and set to sensible default values.
 *
 * @param P pointer to Problem structure
 * @param FCS FCS data structure
 * @param mpi_rank mpi rank of this process
 */
void InitMEMDParameters(struct Problem *P, struct FCSData *FCS, int mpi_rank)
{
  /// set to default values
//  struct memd_params *params = FCS->params.MEMD;

//  const fcs_float timestep = 0;
//  const fcs_int local_number_of_particles = FCS->fcs_particles->n_local_particles;
//  const fcs_float box_size = 1.;
  FCSResult fcs_result;

  Error(SomeError, "fcs_memd_setup not implemented in FCS module, yet.");
//  fcs_result = fcs_memd_setup(FCS->fcs_handle,
//                           box_size,
//                           timestep,
//                           local_number_of_particles,
//                           params->mesh_size,
//                           params->lightspeed,
//                           params->temperature,
//                           params->permittivity
//                           );
//  HandleFCSError(fcs_result);
  fcs_result = fcs_require_virial(FCS->fcs_handle, 1);
  HandleFCSError(fcs_result);
}

/** Frees MEMD specific, allocated memory inside FCSData.
 *
 * @param FCS pointer to FCSData
 */
void FreeMEMDParameters(struct FCSData *FCS)
{
  Free(FCS->params.MEMD);
}

/** Solver specific parser method to get to its required parameters.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param filePos Structure with current position in file
 * @param pd structure for parsed data
 */
void parseMEMDParameters(struct FCSData *FCS, FilePosType *filePos, parse_data *pd)
{
  /// allocate
  FCS->params.MEMD = (struct memd_params *) Malloc( 1*sizeof(struct memd_params), "parseMEMDParameters");
//	struct memd_params *params = FCS->params.MEMD;
//	FCSResult fcs_result;

	while (parse_next_ident(filePos, pd) == 1) {
    switch (pd->ids[0]) {
//    case LC_IDENT(mesh_size):
//      params->mesh_size = parse_int(filePos, pd);
//      fcs_result = fcs_memd_set_mesh_size_1D(FCS->fcs_handle, params->mesh_size);
//      HandleFCSError(fcs_result);
//      break;
//    case LC_IDENT(lightspeed):
//      params->lightspeed = parse_double(filePos, pd);
//      fcs_result = fcs_memd_set_speed_of_light(FCS->fcs_handle, params->lightspeed);
//      HandleFCSError(fcs_result);
//      break;
//    case LC_IDENT(temperature):
//      params->temperature = parse_double(filePos, pd);
//      fcs_result = fcs_memd_set_temperature(FCS->fcs_handle, params->temperature);
//      HandleFCSError(fcs_result);
//      break;
//    case LC_IDENT(permittivity):
//      params->permittivity = parse_double(filePos, pd);
//      fcs_result = fcs_memd_set_permittivity(FCS->fcs_handle, params->permittivity);
//      HandleFCSError(fcs_result);
//      break;
      case LC_IDENT(r_cut):
        FCS->r_cut = parse_double(filePos, pd);
        FCS->r_cut_set = 1;
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
      parse_error(ERRPARSE, filePos, "parseMEMDParameters: invalid variable");
      break;
    }
  }
}

/** Print set of parameters if we are process 0.
 *
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param mpi_rank rank of this process
 */
void printMEMDParameters(struct FCSData *FCS, int mpi_rank)
{
	struct memd_params *params = FCS->params.MEMD;
	FCSResult fcs_result;

  // get all params: MEMD yets misses getter
//  fcs_result = fcs_direct_get_cutoff(FCS->fcs_handle, &params->mesh_size);
//  HandleFCSError(fcs_result);
//  fcs_result = fcs_direct_get_cutoff(FCS->fcs_handle, &params->lightspeed);
//  HandleFCSError(fcs_result);
//  fcs_result = fcs_direct_get_cutoff(FCS->fcs_handle, &params->temperature);
//  HandleFCSError(fcs_result);
//  fcs_result = fcs_direct_get_cutoff(FCS->fcs_handle, &params->permittivity);
//  HandleFCSError(fcs_result);
  fcs_result = fcs_get_tolerance(FCS->fcs_handle, &FCS->tolerance_type, &FCS->tolerance);
  HandleFCSError(fcs_result);

  // print params
  if( mpi_rank == 0 ){
    printf( "Active parameters:\n" );
    printf( "------------------\n" );
    printf( "mesh_size           = %15d\n", params->mesh_size );
    printf( "lightspeed          = %.5E\n", params->lightspeed );
    printf( "temperature         = %.5E\n", params->temperature );
    printf( "permittivity        = %.5E\n", params->permittivity );
    printf( "tolerance           = %.5E\n", FCS->tolerance );
    printf( "tolerance_type      = %d\n", FCS->tolerance_type );
  }
}

/** Method to sum up local energies from all process.
 *
 * @param P pointer to Problem structure
 * @param FCS pointer to FCSData which also contains solver specific parameters
 * @param e_sum_local local sum of energy
 * @param mpi_rank rank of this process
 */
void setMEMDTotalEnergy(struct Problem *P, struct FCSData *FCS, double e_sum_local, int mpi_rank)
{
//	struct memd_params *params = FCS->params.MEMD;
  int i;
  P->Coul->E_l = e_sum_local;

  fcs_float v[9];
  fcs_get_virial(FCS->fcs_handle, v);

  v[i] *= P->epsilon0inv;


  printf("  virial: %e %e %e %e %e %e %e %e %e\n",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
  // inexpensive error estimator is tr(virial) (see spme.c)
  P->Coul->E_err = v[0] + v[4] + v[8];
}

