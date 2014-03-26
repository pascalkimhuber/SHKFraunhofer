/***************************************************************************
                          solparallel_data.cpp  -  description
                             -------------------
    begin                : Wed Jun 08 2004
			   Written by Daniel Bloemer
    email                : bloemer@iam.uni-bonn.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   Tremolo-X Project                                                     *
 *   All rights reserved                                                   *
 *   For details please see the LICENSE file                               *
 *                                                                         *
 ***************************************************************************/

#include "solparallel_data.h"
#include "simulationparameterdata.h"
#include "tremologuiconfig.h"

// QT-Includes
#include <QRegExp>
#include <QString>
#include <QStringList>

#include <cmath>

using std::acos;

#ifdef M_PI
#define PI M_PI
#else
# define PI (acos(-1.0))
#endif

SolParallel_Data::SolParallel_Data( simulationParameterData* parent, TremoloGUIConfig* ConfigIn )
	: r_cut(size_count), 
	  interpolationDegree(size_count), 
	  fcs_tolerance(size_count), 
	  fcs_tolerance_type(size_count)
	  fcs_dipole_correction(size_count), 
	  fcs_max_iterations(size_count), 
	  
{
	// Set parent and gui configuration.
	myParent = parent;
	myConfig = ConfigIn;

	// Set keyword lists, poissonsolver and LB function.
	myKeywordList = new QStringList( QString( "lcs,coulomb,parallelization" ).split(",") );
	poissonsolverList = new QStringList( QString( "FFT,multigrid" ).split(",") );
	poissonsolverFileList = new QStringList( QString( "fft,multigrid" ).split(",") );

	poissonsolver = new QString(poissonsolverList->first());
	
	loadBalancingFunctionList = new QStringList( QString( "OFF,LINEAR,QUADRATIC" ).split(",") );
	loadBalancingFunction = new QString( loadBalancingFunctionList->first() );

	// Set connections to parent simulationParameterData object.
	connect( myParent, SIGNAL( useDynamics( bool ) ), this, SLOT( change_Integration( bool ) ) );
	connect( myParent, SIGNAL( optionHasChanged( whoIsWho ) ), this, SLOT( setAnyOption( whoIsWho ) ) );

	// Set default values. 
	clear();
}

SolParallel_Data::~SolParallel_Data()
{

	// Free memory. 
	delete myKeywordList;
	delete poissonsolverList;
	delete poissonsolverFileList;
	delete poissonsolver;
	delete loadBalancingFunctionList;
	delete loadBalancingFunction;
}

QStringList SolParallel_Data::getKeyWords() const
{
	return *myKeywordList;
}

QString SolParallel_Data::toString() const
{
    // enum longrange_algorithms { OFF = 0, N_2, N_2Spline, SPME, FMM, fcs_direct, fcs_Ewald, fcs_FMM, fcs_PEPC, fcs_PP3MG, fcs_VMG, size_count };

	// Initialize output string. 
	QString t( "" );
	
	// Output selected solver.
	switch(myLongrangeAlgo)
	{
	case OFF: t.append( "No Solver selected\n" );
		break;
	case N_2: t.append( "Solver N^2 selected\n" );
		break;
	case N_2Spline: t.append( "Solver N^2Spline selected\n" );
		break;
	case SPME: t.append( "Solver SPME selected\n" );
		break;
	case FMM: t.append( "Solver FMM selected\n" );
		break;
	case fcs_direct: t.append( "Solver direct (ScaFaCoS) selected\n" );
		break;
	case fcs_Ewald: t.append( "Solver Ewald (ScaFaCoS) selected\n" );
		break;
	case fcs_FMM: t.append( "Solver FMM (ScaFaCoS) selected\n" );
		break;
	case fcs_PEPC: t.append( "Solver PEPC (ScaFaCoS) selected\n" );
		break;
	case fcs_PP3MG: t.append( "Solver PP3MG (ScaFaCoS) selected\n" );
		break;
	case fcs_VMG: t.append( "Solver VMG (ScaFaCoS) selected\n" );
		break;
	}

	// Output selected Poisson solver. 
	t.append( "Poisson solver: " + *poissonsolver + " selected.\n" );

	// Output linked cell r_cut.
	t.append( "Linked Cell r_cut: " + QString::number( linkedCellR_Cut.getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( linkedCellR_Cut ) + "\n" );

	// Output coulomb force constant. 
	t.append( "Permittivity epsilon0inv: " + QString::number( epsilon0inv.getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( epsilon0inv ) + "\n" );

	// Output r_cut of currently selected solver. 
	t.append( "r_cut: " + QString::number( r_cut.at(myLongrangeAlgo).getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( r_cut.at(myLongrangeAlgo) ) + "\n" );

	// Output r_l of currently selected solver. 
	t.append( "r_l: " + QString::number( r_l.getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( r_l) + "\n" );

	// Output splitting coefficient of currently selected solver. 
	t.append( "Splitting Coefficient G: " + QString::number( splittingCoefficientG.getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( splittingCoefficientG ) + "\n" );
	
	// Output MAP of currently selected solver. 
	t.append( "MAP: " + QString::number( MAP.getValue(), 'g' ) + " " + myConfig->trDisplayUnits_SI( MAP ) + "\n" );
	
	// Output cellratio of currently selected solver. 
	t.append( "Cellratio: " + QString::number( isCellratio() ) + "\n" );

	// Output interpolation degree of currently selected solver.
	t.append( "InterpolationDegree: " + QString::number( isInterpolationDegree() ) + "\n" );

	// Output maximum tree level of currently selected solver. 
	t.append( "MaxTreeLevel: " + QString::number( isMaxTreeLevel() ) + "\n" );

	// Output tolerance of currently selected solver. 
	t.append( "Tolerance (ScaFaCoS): " + QString::number( fcs_tolerance.at(myLongrangeAlgo).getValue(), 'g') + " " + myConfig->trDisplayUnits_SI( fcs_tolerance.at(myLongrangeAlgo) ) + "\n" );

	// Output tolerance type of selected solver.
	t.append( "Tolerance type (ScaFaCoS): " + QString::number( isFCS_tolerance_type() ) + "\n" );

	// Output splitting parameter for Ewald method (ScaFaCoS). 
	t.append( "Splitting parameter (Ewald, ScaFaCoS): " + QString::number( fcs_splittingCoefficientAlpha.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI( fcs_splittingCoefficientAlpha) + "\n" );

	// Output K-space cut-off for Ewald method (ScaFaCoS). 
	t.append( "K-space cut-off (Ewald, ScaFaCoS): " + QString::number( isFCS_kmax()) + "\n" );

	// Output maximum K-space cut-off for Ewald method (ScaFaCoS). 
	t.append( "Maximum K-space cut-off (Ewald, ScaFaCoS): " + QString::number (isFCS_maxkmax()) + "\n" );

	// Output number of periodic images used in x-,y- and z-dimension (ScaFaCoS). 
	t.append( "Number of periodic images in x-,y-,z-dimension (ScaFaCoS): " + QString::number( isFCS_periodic_images_x()) + ", " + QString::number( isFCS_periodic_images_y()) + ", " + QString::number( isFCS_periodic_images_z()) + "\n" );

	// Output load balancing routine for FMM (ScaFaCoS).
	t.append( "Load balancing routine for FMM (ScaFaCoS): " + QString::number(isFCS_balanceload()) + "\n");

	// Output dipole correction (FMM and PEPC, ScaFaCoS).
	t.append( "Dipole correction (FMM and PEPC, ScaFaCoS): " + QString::number(isFCS_dipole_correction()) + "\n" );

	// Output maximum tree depth for FMM (ScaFaCoS). 
	t.append( "Maximum tree depth for FMM (ScaFaCoS): " + QString::number(isFCS_maxdepth()) + "\n" );

	// Output potential used for FMM (ScaFaCoS). 
	t.append( "Potential for FMM (ScaFaCoS): " + QString::number(isFCS_potential()) + "\n" );

	// Output radius for CUSP potential in FMM (ScaFaCoS). 
	t.append( "Radius for CUPS potential (FMM, ScaFaCoS): " + QString::number(fcs_radius.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI(fcs_radius) + "\n");

	// Output limit for unrolled functions in FMM (ScaFaCoS). 
	t.append( "Limit for unrolled functions (FMM, ScaFaCoS): " + QString::number(isFCS_unroll_limit()) + "\n");

	// Output debug level for PEPC (ScaFaCoS). 
	t.append( "Debug level (PEPC, ScaFaCoS): " + QString::number(isFCS_debuglevel()) + "\n");

	// Output epsilon for PEPC (ScaFaCoS). 
	t.append("Epsilon (PEPC, ScaFaCoS): " + QString::number(fcs_epsilon.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI(fcs_epsilon) + "\n");

	// Output load balancing switch for PEPC (ScaFaCoS). 
	t.append( "Load balancing switch (PEPC, ScaFaCoS): " + QString::number(isFCS_load_balancing()) + "\n");

	// Output NPM for PEPC (ScaFaCoS). 
	t.append( "NPM (PEPC, ScaFaCoS): " + QString::number(fcs_npm.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI(fcs_npm) + "\n"); 

	// Output number of walk threads per MPI rank for PEPC (ScaFaCoS). 
	t.append("Number of walk threads per MPI (PEPC, ScaFaCoS): " + QString::number(isFCS_num_walk_threads()) + "\n");

	// Output theta for PEPC (ScaFaCoS). 
	t.append("Theta (PEPC, ScaFaCoS): " + QString::number(fcs_theta.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI(fcs_theta) + "\n");

	// Output interpolation degree for PP3MG (ScaFaCoS). 
	t.append("Interpolation degree (PP3MG, ScaFaCoS): " + QString::number(isFCS_degree()) + "\n");

	// Output number of ghost cells for PP3MG (ScaFaCoS). 
	t.append("Number of ghost cells (PP3MG, ScaFaCoS): " + QString::number(isFCS_ghosts()) + "\n");

	// Output number of grid cells in x-,y-,z-dimension for PP3MG (ScaFaCoS).
	t.append("Number of grid cells in x-, y-, z-dimension (PP3MG, ScaFaCoS): " + QString::number(isFCS_gridsize_x()) + ", " + QString::number(isFCS_gridsize_y()) + ", " + QString::number(isFCS_gridsize_z()) + "\n");

	// Output number of iterations for PP3MG and VMG (ScaFaCoS). 
	t.append("Number of iterations (PP3MG and VMG, ScaFaCoS): " + QString::number(isFCS_max_iterations()) + "\n");
	
	// Output multi-grid cycle type specifier for VMG (ScaFaCoS). 
	t.append("Multi-grid cycle type specifier (VMG, ScaFaCoS): " + QString::number(isFCS_cycle_type()) + "\n");

	// Output order of discretization scheme for VMG (ScaFaCoS). 
	t.append("Order of discretization scheme for PDE solver (VMG, ScaFaCoS): " + QString::number(isFCS_discretization_order()) + "\n");

	// Output degree of interpolation polynomial for VMG (ScaFaCoS). 
	t.append("Degree of interpolation polyonmial (VMG, ScaFaCoS): " + QString::number(isFCS_interpolation_order()) + "\n");

	// Output maximum number of multi-grid levels for VMG (ScaFaCoS). 
	t.append("Maximum number of multi-grid levels (VMG, ScaFaCoS): " + QString::number(isFCS_max_level()) + "\n");

	// Output range of smearing of charges over near field cells for VMG (ScaFaCoS). 
	t.append("Smearing range of charges over near field cells (VMG, ScaFaCoS): " + QString::number(isFCS_near_field_cells()) + "\n");

	// Output threshold for residual for VMG (ScaFaCoS). 
	t.append("Residual threshold for solver iteration (VMG, ScaFaCoS): " + QString::number(fcs_precision.getValue(), 'g') + " " + myConfig->trDisplayUnits_SI(fcs_precision) + "\n");

	// Output number of smoothing steps for VMG (ScaFaCoS). 
	t.append("Number of smoothing steps (VMG, ScaFaCoS): " + QString::number(isFCS_smoothing_steps()) + "\n");

	// Output number of processes. 
	if (useSingleProcs) 
	{
		t.append( "Cuts: X: " + QString::number( numberOfProcs[0] ) + " Y: " + QString::number( numberOfProcs[1] ) + " Z: " 
				  + QString::number( numberOfProcs[2] ) + "\tNumber of processors used: " + QString::number( isProcessorUsage() ) + "\n" );
	}
	else 
	{
		t.append( "Processes to Use:\t"  + QString::number( procsAtAll ) + "\n" );
	}
	
	return t;
}

QString SolParallel_Data::toParameterFileString() const
{
	// Initialize empty output string.
	QString t( "" );

	// Output linked cell r_cut. 
	t.append( "lcs:\t" 
			  + writeParamValue( "cellrcut", myConfig->trOut(linkedCellR_Cut) ) + ";\n" );

	// Begin "coulomb" section.
	t.append( "coulomb\t{\n" );

	// Output coulomb force constant. 
	t.append( "\t\tpermittivity:\t" 
			  + writeParamValue( "epsilon0inv", myConfig->trOut(epsilon0inv) ) + ";\n" );

	// Output parameters of N^2 solver.
	t.append( "\t\tn2:\t" 
			  + writeParamState( (myLongrangeAlgo == N_2) ) + ",\t" 
			  + writeParamValue( "r_cut", myConfig->trOut(r_cut.at(N_2)) ) + ",\t" 
			  + writeParamValue( "i_degree", interpolationDegree.at(N_2) ) + ";\n" );

 	// Output parameters of N^2Spline solver. 
	t.append( "\t\tn2spline:\t" 
			  + writeParamState( (myLongrangeAlgo == N_2Spline) ) + ",\t" 
			  + writeParamValue( "r_cut", myConfig->trOut(r_cut.at(N_2Spline)) ) + ",\t" 
			  + writeParamValue( "r_l", myConfig->trOut(r_l.at(N_2Spline)) ) + ",\t" 
			  + writeParamValue( "i_degree", interpolationDegree.at(N_2Spline) ) + ";\n" );

	// Output parameters of TSPME solver.
	t.append( "\t\tspme:\t" 
			  + writeParamState( (myLongrangeAlgo == SPME) ) + ",\t" 
			  + writeParamValue( "ps", (*(poissonsolver.at(SPME))).toLower() ) + ",\t" 
			  + writeParamValue( "r_cut", myConfig->trOut(r_cut.at(SPME)) ) + ",\t" 
			  + writeParamValue( "G", myConfig->trOut(splittingCoefficientG.at(SPME)) ) + ",\t" 
			  + writeParamValue( "cellratio", cellratio.at(SPME) ) + ",\t" 
			  + writeParamValue( "i_degree", interpolationDegree.at(SPME) ) + ";\n" );

	// Output parameters of TFMM solver. 
	t.append( "\t\tfmm:\t" 
			  + writeParamState( (myLongrangeAlgo == FMM) ) + ",\t" 
			  + writeParamValue( "r_cut", myConfig->trOut(r_cut.at(FMM)) ) + ",\t" 
			  + writeParamValue( "i_degree", interpolationDegree.at(FMM) ) + ",\t" 
			  + writeParamValue( "mtl", maxTreeLevel.at(FMM) ) + ",\t" 
			  + writeParamValue( "Map", myConfig->trOut(MAP.at(FMM)) ) + ";\n" );

	// Output parameters of direct summation (ScaFaCoS). 
	t.append( "\t\tfcs_direct:\t"
			  + writeParamState((myLongrangeAlgo == fcs_direct))  + ",\t"
			  + writeParamValue("r_cut", myConfig->trOut(r_cut.at(fcs_direct))) + ",\t"
			  + writeParamValue("periodic_images_x", fcs_periodic_images_x) + ",\t"
			  + writeParamValue("periodic_images_y", fcs_periodic_images_y) + ",\t"
			  + writeParamValue("periodic_images_z", fcs_periodic_images_z) + ";\n"	);

	// Output parameters of Ewald solver (ScaFaCoS). 
	t.append( "\t\tfcs_ewald:\t" 
			  + writeParamState((myLongrangeAlgo == fcs_Ewald) ) + ",\t" 
			  + writeParamValue("r_cut", myConfig->trOut(r_cut.at(fcs_Ewald)) ) + ",\t" 
			  + writeParamValue("tolerance", myConfig->trOut(fcs_tolerance.at(fcs_Ewald))) + ",\t"
			  + writeParamValue("tolerance_type", fcs_tolerance_type.at(fcs_Ewald)) + ",\t" 
			  + writeParamValue("alpha", myConfig->trOut(fcs_splittingCoefficientAlpha)) + ",\t"
			  + writeParamValue("kmax", fcs_kmax) + ",\t"
			  + writeParamValue("maxkmax", fcs_maxkmax) + ";\n" );

	// Output parameters of FMM solver (ScaFaCoS). 
	t.append( "\t\tfcs_fmm:\t"
			  + writeParamState((myLongrangeAlgo == fcs_FMM) ) + ",\t"
			  + writeParamValue("r_cut", myConfig->trOut(r_cut.at(fcs_FMM))) + ",\t"
			  + writeParamValue("tolerance", myConfig->trOut(fcs_tolerance.at(fcs_FMM))) + ",\t"
			  + writeParamValue("tolerance_type", fcs_tolerance_type.at(fcs_FMM)) + ",\t" 
			  + writeParamValue("balanceload", (long) fcs_balanceload) + ",\t"
			  + writeParamValue("dipole_correction", fcs_dipole_correction.at(fcs_FMM)) + ",t"
			  + writeParamValue("maxdepth", (long) fcs_maxdepth) + ",\t"
			  + writeParamValue("potential", fcs_potential) + ",\t"
			  + writeParamValue("radius", myConfig->trOut(fcs_radius)) + ",\t"
			  + writeParamValue("unroll_limit", (long) fcs_unroll_limit) + ";\n" );


	// Close "coulomb" section.
	t.append( "\t\t};\n" );

	// Begin "parallelization" section.
	t.append( "parallelization\t{\n" );

	// Output number of processes.
	t.append( "\t\tproc:\t" 
			  + writeParamValue( "n1", numberOfProcs[0] ) + ",\t" 
			  + writeParamValue( "n2", numberOfProcs[1] ) + ",\t" 
			  + writeParamValue( "n3", numberOfProcs[2] ) + ",\t" 
			  + writeParamValue( "all", procsAtAll ) + ";\n" );

	// Output LB function. 
	t.append( "\t\tloadbal:\t" 
			  + writeParamValue( "wf", (*loadBalancingFunction).toLower() ) + ",\t" 
			  + writeParamValue( "Lb_t", myConfig->trOut(loadbal_t) ) + ",\t" 
			  + writeParamValue( "Lb_delta", myConfig->trOut(loadbal_delta) ) + ",\t" 
			  + writeParamValue( "Lb_step", loadbal_step ) + ";\n" ); 

	// Close "parallelization" section.
	t.append( "\t\t};\n" );

	// Return output string. 
	return t;
}

bool SolParallel_Data::saveValues( QStringList keywordsIn, QStringList identifierListIn, ParseReturnList returnValueListIn )
{
	// Control boolean used to detect errors while parsing.
	bool ok = false;

	// Helper variables used to store temporary values. 
	double d = 0.0;
	int i=0,pos=0;
	QString qs;

	// List of strings containing the names of all longrange solvers. 
	QStringList ensemblelist( QString( "n2,n2spline,spme,fmm,fcs_direct,fcs_ewald,fcs_fmm,fcs_pepc, fcs_pp3mg_pmg,fcs_vmg" ).split(",") );

	// Check if first input keyword is contained in myKeywordList.
	if (myKeywordList->contains(keywordsIn[0].toLower()))
	{
		// Get index of first input keyword in myKeywordList.
     	switch (myKeywordList->indexOf( keywordsIn[0].toLower() ))
		{
		// "LCS" block.
		case 0:
			// Check for "cellrcut" identifier.
			if (identifierListIn[0].toLower() == "cellrcut")
			{
				// Intialize control boolean.
				ok = false;
				
				// Get value of "cellrcut" identifier.
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If a valid value was returned save it to this SolParallel_Data object, else return with error message.
				if (ok) linkedCellR_Cut=myConfig->trIn(linkedCellR_Cut, d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			}
			break;
		// "Coulomb" block. 
 		case 1:
			// "permittivity" section.
			if (keywordsIn[1].toLower() == "permittivity" &&  identifierListIn[0].toLower() == "epsilon0inv")
			{
				// Initialize control boolean.
				ok = false;
				
				// Get value of "epsilon0inv" identifier.
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message. 
				if (ok) epsilon0inv=myConfig->trIn(epsilon0inv, d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Sections of longrange solver: Check for solver and check if its state is "on". 
			else if (ensemblelist.contains( keywordsIn[1].toLower() ) && getParamState(returnValueListIn[0]))
			{
				// + 1 is nescessary, since the enum contains "off". 
				setLongrangeAlgo( ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1 );
			} 
			// Save r_cut identifier.
			else if ( identifierListIn[0].toLower() == "r_cut" ) 
			{
				// Initialize control boolean.
				ok = false;
				
				// Get value of r_cut. 
				d = getParamDoubleValue( returnValueListIn[0], &ok );
				
				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) r_cut.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = myConfig->trIn(r_cut.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1), d ); 
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save r_l identifier.
			else if ( identifierListIn[0].toLower() == "r_l" ) 
			{
				// Initialize control boolean.
				ok = false;

				// Get value of r_l.
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) r_l.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = myConfig->trIn(r_l.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1), d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save splitting coefficient G identifier.
			else if ( identifierListIn[0].toLower() == "g" ) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value of splitting coefficient. 
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) splittingCoefficientG.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) =myConfig->trIn(splittingCoefficientG.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1), d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save cellratio identifier.
			else if ( identifierListIn[0].toLower() == "cellratio" ) 
			{
				// Initialize control boolean.
				ok = false;

				// Get value of cellratio. 
				i = getParamIntValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message. 
				if (ok) cellratio.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save interpolation degree.
			else if ( identifierListIn[0].toLower() == "i_degree" ) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value of interpolation degree. 
				i = getParamIntValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) interpolationDegree.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save poisson solver. 
			else if ( identifierListIn[0].toLower() == "ps" ) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get QString for the poisson solver. 
				qs = getParamStringValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) { if (poissonsolverFileList->contains( qs.toLower() )) *(poissonsolver.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1)) = (*poissonsolverList)[poissonsolverFileList->indexOf( qs.toLower() )]; }
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save maximum tree level. 
			else if ( identifierListIn[0].toLower() == "mtl" ) 
			{
				// Initialize control boolean. 
				ok = false;
				
				// Get value for maximum tree level. 
				i = getParamIntValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) maxTreeLevel.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save MAP. 
			else if ( identifierListIn[0].toLower() == "map" ) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value for MAP. 
				d = getParamDoubleValue( returnValueListIn[0], &ok );

                // If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) MAP.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1) = myConfig->trIn(MAP.at(ensemblelist.indexOf( keywordsIn[1].toLower() ) + 1), d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Do nothing if no longrange solver is "on". 
			else if (ensemblelist.contains( keywordsIn[1].toLower() ) && !getParamState(returnValueListIn[0])) 
			{
				// do nothing, but don't produce an error.
			} 
			// Return error whenever an invalid parameter is parsed. 
			else  
			{ 
				// Returned with error message. 
				myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "" ); 
				return false; 
			}
			break;
		// "Parallelization" block. 
		case 2:		
			// Save number of processes in every dimension. 
			if ((keywordsIn[1].toLower() == "proc") && (identifierListIn[0].toLower().startsWith( "n" ))) 
			{
				// Get position (n1, n2 or n3). 
				pos = identifierListIn[0].toLower().remove("n").toInt();

				// Initialize control boolean. 
				ok = false;

				// Get value for the number of processes. 
				i = getParamIntValue( returnValueListIn[0], &ok );

                // If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if ((pos>0) && (pos<4) && (ok) && (i>0)) numberOfProcs[pos-1] = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save whole number of processes.
			else if ((keywordsIn[1].toLower() == "proc") && (identifierListIn[0].toLower() == "all" )) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get number of processes. 
				i = getParamIntValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if ((ok) && (i>0)) procsAtAll = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save LB function. 
			else if ((keywordsIn[1].toLower() == "loadbal") && (identifierListIn[0].toLower() == "wf")) 
			{
				// Initialize control boolean.
				ok = false;

				// Get QString for the LB functions. 
				qs = getParamStringValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) setLoadBalancingFunction( qs.upper() );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save LB starttime.
			else if ((keywordsIn[1].toLower() == "loadbal") && (identifierListIn[0].toLower() == "lb_t")) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value for LB starttime. 
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) loadbal_t=myConfig->trIn(loadbal_t, d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save LB timedelta
			else if ((keywordsIn[1].toLower() == "loadbal") && (identifierListIn[0].toLower() == "lb_delta")) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value for LB timedelta.
				d = getParamDoubleValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) loadbal_delta=myConfig->trIn(loadbal_delta,d );
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Save LB step. 
			else if ((keywordsIn[1].toLower() == "loadbal") && (identifierListIn[0].toLower() == "lb_step")) 
			{
				// Initialize control boolean. 
				ok = false;

				// Get value for LB step. 
				i = getParamIntValue( returnValueListIn[0], &ok );

				// If returned value was valid save it to this SolParallel_Data object, else return with error message.
				if (ok) loadbal_step = i;
				else { myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "invalid format" );  return false; }
			} 
			// Return error whenever an invalid parameter is parsed. 
			else  
			{ 
				// Return with error message.
				myConfig->addLogEntryInvalidParm( TremoloGUIConfig::error, returnValueListIn[0], "" );  
				return false; 
			}
			break;
		}
	}

	// Return true if no error occured. 
	return true;
}

double SolParallel_Data::isLinkedCellR_Cut() 
{ return myConfig->trOut( linkedCellR_Cut ); }

QString SolParallel_Data::isLinkedCellR_CutUnit() 
{ return myConfig->trDisplayUnits( linkedCellR_Cut ); }

double SolParallel_Data::isEpsilon0inv()
{ return myConfig->trOut( epsilon0inv ); }

QString SolParallel_Data::isEpsilon0invUnit()
{ return myConfig->trDisplayUnits( epsilon0inv ); }

double SolParallel_Data::isR_Cut() 
{ return myConfig->trOut( r_cut.at(myLongrangeAlgo) ); }

QString SolParallel_Data::isR_CutUnit() 
{ return myConfig->trDisplayUnits( r_cut.at(myLongrangeAlgo) ); }

double SolParallel_Data::isR_L() 
{ return myConfig->trOut( r_l ); }

QString SolParallel_Data::isR_LUnit() 
{ return myConfig->trDisplayUnits( r_l ); }

double SolParallel_Data::isSplittingCoefficientG() 
{ return myConfig->trOut( splittingCoefficientG ); }

QString SolParallel_Data::isSplittingCoefficientGUnit() 
{ return myConfig->trDisplayUnits( splittingCoefficientG ); }

double SolParallel_Data::isMAP() 
{ return myConfig->trOut( MAP ); }

QString SolParallel_Data::isMAPUnit() 
{ return myConfig->trDisplayUnits( MAP ); }

double SolParallel_Data::isFCS_tolerance()
{ return myConfig->trOut( fcs_tolerance.at(myLongrangeAlgo) ); }

QString SolParallel_Data::isFCS_toleranceUnit()
{ return myConfig->trDisplayUnits( fcs_tolerance.at(myLongrangeAlgo) ); }

double SolParallel_Data::isFCS_splittingCoefficientAlpha()
{ return myConfig->trOut( fcs_splittingCoefficientAlpha); }

double SolParallel_Data::isFCS_radius()
{ return myConfig->trOut( fcs_radius); }

QString SolParallel_Data::isFCS_radiusUnit()
{ return myConfig->trDisplayUnits( fcs_radius); }

double SolParallel_Data::isFCS_epsilon()
{ return myConfig->trOut( fcs_epsilon); }

QString SolParallel_Data::isFCS_epsilonUnit()
{ return myConfig->trDisplayUnits( fcs_epsilon); }

double SolParallel_Data::isFCS_theta()
{ return myConfig->trOut( fcs_theta); }

QString SolParallel_Data::isFCS_thetaUnit()
{ return myConfig->trDisplayUnits( fcs_theta); }

double SolParallel_Data::isFCS_precision()
{ return myConfig->trOut( fcs_precision); }

QString SolParallel_Data::isFCS_precisionUnit()
{ return myConfig->trDisplayUnits( fcs_precision); }

long SolParallel_Data::isProcessorUsage() const
{ return ( numberOfProcs[0] * numberOfProcs[1] * numberOfProcs[2] ); }

double SolParallel_Data::isLoadbal_t() 
{ return myConfig->trOut( loadbal_t ); }

QString SolParallel_Data::isLoadbal_tUnit() 
{ return myConfig->trDisplayUnits( loadbal_t ); }

double SolParallel_Data::isLoadbal_delta() 
{ return myConfig->trOut( loadbal_delta ); }

QString SolParallel_Data::isLoadbal_deltaUnit() 
{ return myConfig->trDisplayUnits( loadbal_delta ); }

void SolParallel_Data::clear()
{
	// electric constant epsilon0 in F/m: http://physics.nist.gov/cgi-bin/cuu/Value?ep0
	const double epsilon0 = 8.854187817e-12;	// A^2 s^4 / kg m^3

	// Set default values for linked cell r_cut and coulomb force constant.
	linkedCellR_Cut.set(2.5e-10);
	linkedCellR_Cut.set(1,0,0,0,0);
	epsilon0inv.set(1. / (4. * PI * epsilon0));
	epsilon0inv.set(3,-4,1,-2,0);

	// Disable all longrange solvers. 
	setLongrangeAlgo( OFF );

	// Set default values for the parameters of ALL  longrange algorithms.	
	for (int solver = 0; solver < size_count; solver++) {
		
		// Set default poissonsolver for all longrangre solvers. 
		*poissonsolver.at(solver) = poissonsolverList->first();

		// Set default values (and unit exponents) for r_cut, r_l, splittingCoefficientG, MAP, cellratio, interpolationDegree and maxTreeLevel. 
		r_cut.at(solver).set(2.5e-10);
		r_cut.at(solver).set(1, 0, 0, 0, 0);
		
		r_l.at(solver).set(2.2e-10);
		r_l.at(solver).set(1, 0, 0, 0, 0);

		splittingCoefficientG.at(solver).set(0.35e10);
		splittingCoefficientG.at(solver).set(-1, 0, 0, 0, 0);

		MAP.at(solver).set(0.577);
		MAP.at(solver).set(0, 0, 0, 0, 0);

		cellratio.at(solver) = 6;
		interpolationDegree.at(solver) = 5;
		maxTreeLevel.at(solver) = 6;
	}

	// Set default values for parallelization parameters. 
	numberOfProcs[0] = 1;
	numberOfProcs[1] = 1;
	numberOfProcs[2] = 1;
	procsAtAll = 1;
	useSingleProcs = true;
	
	*loadBalancingFunction = loadBalancingFunctionList->first();
	emit enableLoadBalValues_dyn( false );
	emit enableLoadBalValues_opt( false );
	
	loadbal_t.set(0.0);
	loadbal_t.set(0,1,0,0,0);
	loadbal_delta.set(1.0);
	loadbal_delta.set(0,1,0,0,0);
	loadbal_step = 1;
	
	// Emit hasChanged() signal. 
	emit hasChanged();
}

void SolParallel_Data::setLongrangeAlgo( int algo )
{
	switch(algo)
	{
	case OFF: myLongrangeAlgo = OFF;
		emit enablepoisson_solver( false );
		emit enableR_Cut( false );
		emit enableR_L( false );
		emit enableSplittingCoefficientG( false );
		emit enableMAP( false );
		emit enableMaxTreeLevel( false );
		emit enableInterpolationDegree( false );
		emit enableCellratio( false );
		useSingleProcs = true;
		break;
	case N_2: myLongrangeAlgo = N_2;
		emit enablepoisson_solver( false );
		emit enableR_Cut( true );
		emit enableR_L( false );
		emit enableSplittingCoefficientG( false );
		emit enableMAP( false );
		emit enableMaxTreeLevel( false );
		emit enableInterpolationDegree( true );
		emit enableCellratio( false );
		useSingleProcs = true;
		break;
	case N_2Spline: myLongrangeAlgo = N_2Spline;
		emit enablepoisson_solver( false );
		emit enableR_Cut( true );
		emit enableR_L( true );
		emit enableSplittingCoefficientG( false );
		emit enableMAP( false );
		emit enableMaxTreeLevel( false );
		emit enableInterpolationDegree( true );
		emit enableCellratio( false );
		useSingleProcs = true;
		break;
	case Ewald: myLongrangeAlgo = fcs_Ewald;
		emit enablepoisson_solver( false );
		emit enableR_Cut( true );
		emit enableR_L( false );
		emit enableSplittingCoefficientG( true );
		emit enableMAP( false );
		emit enableMaxTreeLevel( false );
		emit enableInterpolationDegree( true );
		emit enableCellratio( true );
		useSingleProcs = true;
		break;
	case SPME: myLongrangeAlgo = SPME;
		emit enablepoisson_solver( true );
		emit enableR_Cut( true );
		emit enableR_L( false );
		emit enableSplittingCoefficientG( true );
		emit enableMAP( false );
		emit enableMaxTreeLevel( false );
		emit enableInterpolationDegree( true );
		emit enableCellratio( true );
		useSingleProcs = true;
		break;
	case FMM: myLongrangeAlgo = FMM;
		emit enablepoisson_solver( false );
		emit enableR_Cut( true );
		emit enableR_L( false );
		emit enableSplittingCoefficientG( false );
		emit enableMAP( true );
		emit enableMaxTreeLevel( true );
		emit enableInterpolationDegree( true );
		emit enableCellratio( false );
		useSingleProcs = false;
		break;
	}
	emit hasChanged();
}

void SolParallel_Data::setLinkedCellR_Cut( double valueIn )
{
	linkedCellR_Cut = myConfig->trIn( linkedCellR_Cut, valueIn );
	emit hasChanged();
}

void SolParallel_Data::setEpsilon0inv( double valueIn )
{
	epsilon0inv = myConfig->trIn( epsilon0inv, valueIn );
	emit hasChanged();
}

void SolParallel_Data::setR_Cut( double valueIn )
{
	r_cut.at(myLongrangeAlgo) = myConfig->trIn( r_cut.at(myLongrangeAlgo), valueIn );
	
	if (! (r_l.at(myLongrangeAlgo) < r_cut.at(myLongrangeAlgo))) 
	{
     	r_l.at(myLongrangeAlgo) = r_cut.at(myLongrangeAlgo) * 0.9;
     	myConfig->addLogEntry( TremoloGUIConfig::warning, "Reducing r_l to be smaller than r_cut." );
	}
	emit hasChanged();
}

void SolParallel_Data::setR_L( double valueIn )
{
	StandardisedDouble temp = myConfig->trIn( r_l.at(myLongrangeAlgo), valueIn );

	if (temp < r_cut.at(myLongrangeAlgo)) 
		r_l.at(myLongrangeAlgo) = temp; 
	else 
		myConfig->addLogEntry( TremoloGUIConfig::warning, "Value for r_l not accepted, since it was bigger than r_cut." );
	emit hasChanged();
}

void SolParallel_Data::setSplittingCoefficientG( double valueIn )
{
	splittingCoefficientG.at(myLongrangeAlgo) = myConfig->trIn( splittingCoefficientG.at(myLongrangeAlgo), valueIn );
	emit hasChanged();
}

void SolParallel_Data::setMAP( double valueIn )
{
	MAP.at(myLongrangeAlgo) = myConfig->trIn( MAP.at(myLongrangeAlgo), valueIn );
	emit hasChanged();
}

void SolParallel_Data::setCellratio( int valueIn )
{
	cellratio.at(myLongrangeAlgo) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setInterpolationDegree( int valueIn )
{
	interpolationDegree.at(myLongrangeAlgo) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setMaxTreeLevel( int valueIn )
{
	maxTreeLevel.at(myLongrangeAlgo) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setPoisson_solver( const QString &valueIn )
{
	*(poissonsolver.at(myLongrangeAlgo)) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_tolerance(double valueIn)
{ 
	fcs_tolerance.at(myLongrangeAlgo) = myConfig->trIn( fcs_tolerance.at(myLongrangeAlgo), valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_tolerance_type(int valueIn)
{
	fcs_tolerance_type.at(myLongrangeAlgo) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_splittingCoefficientAlpha(double valueIn)
{
	fcs_splittingCoefficientAlpha.at(myLongrangeAlgo) = myConfig->trIn( fcs_splittingCoefficientAlpha.at(myLongrangeAlgo), valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_kmax(int valueIn)
{
	fcs_kmax = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_maxkmax(int valueIn)
{
	fcs_maxkmax = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_periodic_images_x(int valueIn)
{
	fcs_periodic_images_x = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_periodic_images_y(int valueIn)
{
	fcs_periodic_images_y = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_periodic_images_z(int valueIn)
{
	fcs_periodic_images_z = valueIn;
	emit hasChanged();
}
	
void SolParallel_Data::setFCS_balanceload(int valueIn)
{
	fcs_balanceload = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_dipole_correction(int valueIn)
{
	fcs_dipole_correction = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_maxdepth(long long valueIn)
{
	fcs_maxdepth = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_potential(int valueIn)
{
	fcs_potential = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_radius(double valueIn)
{
	fcs_radius = myConfig->trIn(fcs_radius, valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_unroll_limit(long long valueIn)
{
	fcs_unroll_limit = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_debuglevel(int valueIn)
{
	fcs_debuglevel = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_epsilon(double  valueIn)
{
	fcs_epsilon = myConfig->trIn(fcs_epsilon, valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_load_balancing(int valueIn)
{
	fcs_load_balancing = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_npm(double valueIn)
{
	fcs_npm = myConfig->trIn(fcs_npm, valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_num_walk_threads(int valueIn)
{
	fcs_num_walk_threads = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_theta(double valueIn)
{
	fcs_theta = myConfig->trIn(fcs_theta, valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_degree(int valueIn)
{
	fcs_degree = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_ghosts(int valueIn)
{
	fcs_ghosts = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_gridsize_x(int valueIn)
{
	fcs_gridsize_x = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_gridsize_y(int valueIn)
{
	fcs_gridsize_y = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_gridsize_z(int valueIn)
{
	fcs_gridsize_z = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_max_iterations(int valueIn)
{
	fcs_max_iterations.at(myLongrangeAlgo) = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_cycle_type(int valueIn)
{
	fcs_cycle_type = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_discretization_order(int valueIn)
{
	fcs_discretization_order = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_interpolation_order(int valueIn)
{
	fcs_interpolation_order = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_max_level(int valueIn)
{
	fcs_max_level = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_near_field_cells(int valueIn)
{
	fcs_near_field_cells = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setFCS_precision(double valueIn)
{
	fcs_precision = myConfig->trIn(fcs_precision, valueIn);
	emit hasChanged();
}

void SolParallel_Data::setFCS_smoothing_steps(int valueIn)
{
	fcs_smoothing_steps = valueIn;
	emit hasChanged();
}

void SolParallel_Data::setCutsX( int cutsIn)
{
	if (cutsIn > 0) numberOfProcs[0] = cutsIn;
	emit hasChanged();
}

void SolParallel_Data::setCutsY( int cutsIn)
{
	if (cutsIn > 0) numberOfProcs[1] = cutsIn;
	emit hasChanged();
}

void SolParallel_Data::setCutsZ( int cutsIn)
{
	if (cutsIn > 0) numberOfProcs[2] = cutsIn;
	emit hasChanged();
}
     
void SolParallel_Data::setProcsAtAll( int procsIn )
{
	if (procsIn > 0) procsAtAll = procsIn;
	emit hasChanged();
}
     
void SolParallel_Data::setLoadBalancingFunction( const QString &valueIn )
{
	*loadBalancingFunction = valueIn;
	if (loadBalancingFunctionList->indexOf( *loadBalancingFunction ) > 0) 
		if (myParent->isIntegrationScheme() == simulationParameterData::dynamics) 
		{
			emit enableLoadBalValues_dyn( true );
			emit enableLoadBalValues_opt( false );
		} else {
			emit enableLoadBalValues_dyn( false );
			emit enableLoadBalValues_opt( true );
		} else {
		emit enableLoadBalValues_dyn( false );
		emit enableLoadBalValues_opt( false );
	}
	emit hasChanged();
}

void SolParallel_Data::setLoadbal_t( double valueIn )
{
	loadbal_t = myConfig->trIn( loadbal_t, valueIn );
	emit hasChanged();
}

void SolParallel_Data::setLoadbal_delta( double valueIn )
{
	loadbal_delta = myConfig->trIn( loadbal_delta, valueIn );
	emit hasChanged();
}

void SolParallel_Data::setLoadbal_step( int valueIn)
{
	if (valueIn > 0) loadbal_step = valueIn;
	emit hasChanged();
}

void SolParallel_Data::change_Integration( bool dynamo )
{
	if (loadBalancingFunctionList->indexOf( *loadBalancingFunction ) > 0) 
	{
     	if (dynamo) {
			emit enableLoadBalValues_dyn( true );
			emit enableLoadBalValues_opt( false );
		} else {
			emit enableLoadBalValues_dyn( false );
			emit enableLoadBalValues_opt( true );
		}
	}
}

void SolParallel_Data::setAnyOption( whoIsWho whatShouldISet )
{
	switch (whatShouldISet)
	{
	case simulationParameterdDataTypes::cAll:
		emit hasChanged();
		break;
	}

}

#include "solparallel_data.moc"
