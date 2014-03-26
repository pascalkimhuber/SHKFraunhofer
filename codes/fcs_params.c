/*
 * Project: MoleCuilder
 * Description: creates and alters molecular systems
 * Copyright (C)  2012 University of Bonn. All rights reserved.
 * Please see the LICENSE file or "Copyright notice" in builder.cpp for details.
 */

/*
 * fcs_params.c
 *
 *  Created on: Feb 1, 2012
 *      Author: heber
 */

#include "fcs_params.h"

#include <fcs.h>

#include "parser.h"

#include "direct_params.h"
#include "ewald_params.h"
#include "fmm_params.h"
#include "pepc_params.h"
#include "pp3mg_pmg_params.h"
#include "vmg_params.h"

#ifdef USE_GPL
#include "p2nfft_params.h"
#include "p3m_params.h"
#include "memd_params.h"
#include "mmm_params.h"
#endif

/* order must match enum FCSMethod ! */

#ifdef USE_GPL
const FCSInit_array FCSInit[] = {
	{ LC_KEY(fcs_direct), "DIRECT", "direct", &InitDIRECTParameters, &parseDIRECTParameters, &printDIRECTParameters, &setDIRECTTotalEnergy, &FreeDIRECTParameters },
	{ LC_KEY(fcs_ewald), "EWALD", "ewald", &InitEWALDParameters, &parseEWALDParameters, &printEWALDParameters, &setEWALDTotalEnergy, &FreeEWALDParameters },
	{ LC_KEY(fcs_fmm), "FMM", "fmm", &InitFMMParameters, &parseFMMParameters, &printFMMParameters, &setFMMTotalEnergy, &FreeFMMParameters },
	{ LC_KEY(fcs_memd), "MEMD", "memd", &InitMEMDParameters, &parseMEMDParameters, &printMEMDParameters, &setMEMDTotalEnergy, &FreeMEMDParameters },
	{ LC_KEY(fcs_mmm), "MMM", "mmm", &InitMMMParameters, &parseMMMParameters, &printMMMParameters, &setMMMTotalEnergy, &FreeMMMParameters },
	{ LC_KEY(fcs_pepc), "PEPC", "pepc", &InitPEPCParameters, &parsePEPCParameters, &printPEPCParameters, &setPEPCTotalEnergy, &FreePEPCParameters },
	{ LC_KEY(fcs_p2nfft), "P2NFFT", "p2nfft", &InitP2NFFTParameters, &parseP2NFFTParameters, &printP2NFFTParameters, &setP2NFFTTotalEnergy, &FreeP2NFFTParameters },
	{ LC_KEY(fcs_p3m), "P3M", "p3m", &InitP3MParameters, &parseP3MParameters, &printP3MParameters, &setP3MTotalEnergy, &FreeP3MParameters },
	{ LC_KEY(fcs_pp3mg_pmg), "PP3MG_PMG", "pp3mg", &InitPP3MG_PMGParameters, &parsePP3MG_PMGParameters, &printPP3MG_PMGParameters, &setPP3MG_PMGTotalEnergy, &FreePP3MG_PMGParameters },
  { LC_KEY(fcs_vmg), "VMG", "vmg", &InitVMGParameters, &parseVMGParameters, &printVMGParameters, &setVMGTotalEnergy, &FreeVMGParameters }
};
#else
const FCSInit_array FCSInit[] = {
    { LC_KEY(fcs_direct), "DIRECT", "direct", &InitDIRECTParameters, &parseDIRECTParameters, &printDIRECTParameters, &setDIRECTTotalEnergy, &FreeDIRECTParameters },
    { LC_KEY(fcs_ewald), "EWALD", "ewald", &InitEWALDParameters, &parseEWALDParameters, &printEWALDParameters, &setEWALDTotalEnergy, &FreeEWALDParameters },
    { LC_KEY(fcs_fmm), "FMM", "fmm", &InitFMMParameters, &parseFMMParameters, &printFMMParameters, &setFMMTotalEnergy, &FreeFMMParameters },
    { LC_KEY(fcs_memd), "MEMD", "memd", NULL, NULL, NULL, NULL, NULL },
    { LC_KEY(fcs_mmm), "MMM", "mmm", NULL, NULL, NULL, NULL, NULL },
    { LC_KEY(fcs_pepc), "PEPC", "pepc", &InitPEPCParameters, &parsePEPCParameters, &printPEPCParameters, &setPEPCTotalEnergy, &FreePEPCParameters },
    { LC_KEY(fcs_p2nfft), "P2NFFT", "p2nfft", NULL, NULL, NULL, NULL, NULL },
    { LC_KEY(fcs_p3m), "P3M", "p3m", NULL, NULL, NULL, NULL, NULL},
    { LC_KEY(fcs_pp3mg_pmg), "PP3MG_PMG", "pp3mg", &InitPP3MG_PMGParameters, &parsePP3MG_PMGParameters, &printPP3MG_PMGParameters, &setPP3MG_PMGTotalEnergy, &FreePP3MG_PMGParameters },
  { LC_KEY(fcs_vmg), "VMG", "vmg", &InitVMGParameters, &parseVMGParameters, &printVMGParameters, &setVMGTotalEnergy, &FreeVMGParameters }
};
#endif

/*@constant int N_FCS;@*/
#define N_FCS ((int)(sizeof FCSInit/sizeof FCSbInit[0]))

