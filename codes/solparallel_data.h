/***************************************************************************
                          solparallel_data.h  -  DataStorage for solver and parallelization
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

#ifndef HAVE_SOLPARALLEL
#define HAVE_SOLPARALLEL

// Standard includes.
#include <vector>


#include "tgcommon.h"
#include "standardiseddouble.h"
#include "parameter_template.h"

// Qt includes.
#include <QObject>

class simulationParameterData;
class TremoloGUIConfig;
class QString;
class QStringList;

/**
 * @brief This class stores the data about solver and parallelization in the Simulation parameters section. 
 */
class SolParallel_Data : public Parameter_template, public TremoloGUI_CommonTypes 
{
     Q_OBJECT

public:
	/**
	 * @brief Constructor.
	 *        The constructor initializes all parameter vectors (e.g. r_cut, ...) such that they have all the same length equal to the number of longrange solvers. 
	 *
	 * @param parent Parent simulationParameterData object.
	 * @param configIn Config data for this object. 
	 */
	SolParallel_Data( simulationParameterData* parent, TremoloGUIConfig* configIn );

	/// Destructor. 
	~SolParallel_Data();
     
	/**
	 * @brief Enumeration of all longrange solvers.
	 *        Note that the last element in the enum is only used to determine the number of longrange solvers!
	 */
	enum longrange_algorithms { OFF = 0, N_2, N_2Spline, SPME, FMM, fcs_direct, fcs_Ewald, fcs_FMM, fcs_PEPC, fcs_PP3MG, fcs_VMG, fcs_P2NFFT, size_count };
     
	/**
	 * @brief Returns the index of the longrange solver that is currently active in this SolParallel_Data object.
	 *
	 * @return Index of current longrange solver.
	 */
	longrange_algorithms isLongrangeAlgo() { return myLongrangeAlgo; }

	/**
	 * @brief Returns a list of keywords used by this SolParallel_Data object for parsing the parameters file. 
	 *
	 * @return QStringList containing the keywords used in the parameters file. 
	 */
	QStringList getKeyWords() const;

	/**
	 * @brief Returns debug output of all data contained in this SolParallel_Data object.
	 *
	 * @note Only the data corresponding to the currently selected longrange solver is output. 
	 * @note For parameters that are not used by the currently selected solver (e.g. r_l for N^2  solver) only default values are returned!
	 *
	 * @return QString containing debug output. 
	 */
	QString toString() const;
	 
	/**
	 * @brief Writes the parameters of this SolParallel_Data object to the parameters file. 
	 *
	 * @return QString containing the parameters of the @c coulomb and the @c parallelization section for the parameters file.
	 */
	QString toParameterFileString() const;
	
	/**
	 * @brief Saves the value of the "Solver and Parallelization" section of the parameters file to this SolParallel_Data object.
	 *        The data has been parsed in a previous step and is now passed to this method.
	 *
	 * @param keywordsIn A list of keywords containing the keywords of the "Solver and Parallelization section of the parameters file".
	 * @param identifierListIn A list of identifiers (corresponding to the keywordsIn parameter).
	 * @param returnValueListIn A list containing the values of the corresponding identifiers. 
	 *
	 * @return False if any error occured, true otherwise. 
	 */
	bool saveValues( QStringList keywordsIn, QStringList identifierListIn, ParseReturnList returnValueListIn );

	/**
	 * @brief Returns scaled value of linked cell r_cut. 
	 *
	 * @return Scaled value of linked cell r_cut. 
	 */
	double isLinkedCellR_Cut();

	/**
	 * @brief Returns the units for linked cell r_cut. 
	 *
	 * @return String containing the units for the linked cell r_cut value. 
	 */
	QString isLinkedCellR_CutUnit();

	/**
	 * @brief Returns the scaled value of the coulomb force constant @c epsilon0inv = \f$\frac{1}{4\pi\epsilon_0} \cdot \frac{e^2}{[length]}\f$. 
	 *
	 * @return Scaled value of @c epsilon0inv.  
	 */
	double isEpsilon0inv();

	/**
	 * @brief Returns the units for the coulomb force constant. 
	 *
	 * @return String containing the units for the coulomb force constant.
	 */
	QString isEpsilon0invUnit();

	/**
	 * @brief Returns the scaled r_cut value of the currently selected longrange solver. 
	 *
	 * @return Scaled value of the current solver's r_cut parameter.
	 */
	double isR_Cut();

	/**
	 * @brief Returns the units of the r_cut parameter of the currently selected longrange solver. 
	 *
	 * @return String containing the units for the current solver's r_cut parameter. 
	 */ 
	QString isR_CutUnit();

	/**
	 * @brief Returns the scaled r_l value of the currently selected longrange solver. 
	 * 
	 * @return Scaled value of the current solver's r_l parameter. 
	 */
	double isR_L();

	/**
	 * @brief Returns the units for the r_l parameter of the currently selected longrange solver. 
	 *
	 * @return String containing the units for the current solver's r_l parameter. 
	 */
	QString isR_LUnit();

	/**
	 * @brief Returns the scaled value of the splitting coefficient G of the currently selected longrange solver. 
	 * 
	 * @return Scaled value of the current solver's splitting coefficient. 
	 */	
	double isSplittingCoefficientG();

	/**
	 * @brief Returns the units for the splitting coefficient of the currently selected longrange solver. 
	 *
	 * @return String containing the units for the current solver's splitting coefficient.
	 */
	QString isSplittingCoefficientGUnit();

	/**
	 * @brief Returns the scaled value of the MAP parameter of the currently selected longrange solver. 
	 * 
	 * @return Scaled value of the current solver's MAP parameter. 
	 */
	double isMAP();

	/**
	 * @brief Returns the units for the MAP parameter of the currently selected longrange solver. 
	 *
	 * @return String containing the units for the current solver's MAP parameter.
	 */
	QString isMAPUnit();

	/// @warning Is not implemented?
	int isGitterOrCell();

	/// @warning Is not implemented?
	QString isGitterOrCellUnit();

	/**
	 * @brief Returns the interpolation degree of the currently selected longrange solver. 
	 *
	 * @return Interpolation degree of current solver.
	 */
	int isInterpolationDegree() const { return interpolationDegree.at(myLongrangeAlgo); }
	/**
	 * @brief Returns the maximum tree level of the currently selected longrange solver. 
	 *
	 * @return Maximum tree level of current multipole method. 
	 */
	int isMaxTreeLevel() const { return maxTreeLevel; }

	/**
	 * @brief Returns the list of available poisson solver for the currently selected longrange algorithm.  
	 *
	 * @return QStringList containing all available poisson solver for the current solver.
	 */
	QStringList* isPoisson_solverList() const { return poissonsolverList; }

	/**
	 * @brief Returns the name of the selected poissonsolver for the current longrange solver. 
	 *
	 * @return QString containing the name of the actual selected poissonsolver.
	 */
	QString isPoisson_solver() const { return *poissonsolver; }

	/**
	 * @brief Returns the cellratio of the currently selected longrange solver. 
	 *
	 * @return Cellratio of current solver.
	 */
	int isCellratio() const { return cellratio; }

	/**
	 * @brief Returns the tolerance (for ScaFaCos solvers only) for the current longrange solver. 
	 *
	 * @return Tolerance (see ScaFaCoS) of the current solver. 
	 */
	double isFCS_tolerance();

	/**
	 * @brief Returns the units for the tolerance (for ScaFaCoS solvers only) for the current longrange solver. 
	 *
	 * @return QString containing the units for the current solver's tolerance. 
	 */
	QString isFCS_toleranceUnit();

	/**
	 * @brief Returns the tolerance type (for ScaFaCoS solvers only) for the current longrange solver.
	 * 
	 * @return Tolerance type (see ScaFaCos) of current solver. 
	 */
	int isFCS_tolerance_type() const { return fcs_tolerance_type.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the splitting coefficient for the Ewald method and P2NFFT (ScaFaCoS). 
	 *
	 * @return Splitting coefficient \f$\alpha\f$. 
	 */
	double isFCS_splittingCoefficientAlpha();

	/**
	 * @brief Returns the units for the splitting coefficient for the Ewald method and P2NFFT (ScaFaCoS). 
	 *
	 * @return QString containing the units for the splitting coefficient \f$\alpha\f$. 
	 */
	QString isFCS_splittingCoefficientAlphaUnit();

	/**
	 * @brief Returns the K-space cut-off for the Ewald method (ScaFaCoS). 
	 *
	 * @return K-space cut-off for the Ewald method. 
	 */
	int isFCS_kmax() const { return fcs_kmax; }

	/**
	 * @brief Returns the maximum K-space cut-off for the Ewald method (ScaFaCoS).
	 *
	 * @return Maximum K-space cut-off for the Ewald method. 
	 */
	int isFCS_maxkmax() const { return fcs_maxkmax; }

	/**
	 * @brief Returns the number of periodic images used in x-dimension for the direct method (ScaFaCoS). 
	 *
	 * @return Number of periodic images used in x-dimension for the direct method. 
	 */
	int isFCS_periodic_images_x() const { return fcs_periodic_images_x; }

	/**
	 * @brief Returns the number of periodic images used in y-dimension for the direct method (ScaFaCoS). 
	 *
	 * @return Number of periodic images used in y-dimension for the direct method. 
	 */
	int isFCS_periodic_images_y() const { return fcs_periodic_images_y; }

	/**
	 * @brief Returns the number of periodic images used in z-dimension for the direct method (ScaFaCoS). 
	 *
	 * @return Number of periodic images used in z-dimension for the direct method. 
	 */
	int isFCS_periodic_images_z() const { return fcs_periodic_images_z; }

	/**
	 * @brief Returns the index of the load balancing routine for FMM (ScaFaCoS).
	 *
	 * @return Index of the load balancing routine for FMM. 
	 */
	int isFCS_balanceload() const { return fcs_balanceload; }

	/**
	 * @brief Returns the dipole correction for FMM and PEPC (ScaFaCoS). 
	 *
	 * @return Dipole correction. 
	 */
	int isFCS_dipole_correction() const { return fcs_dipole_correction.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the maximum tree depth for FMM (ScaFaCoS). 
	 *
	 * @return Maximum tree depth for FMM. 
	 */
	long long isFCS_maxdepth() const { return fcs_maxdepth; }

	/**
	 * @brief Returns the index of the potential used for FMM (ScaFaCoS). 
	 *
	 * @return Index of the potential for FMM. 
	 */
	int isFCS_potential() const { return fcs_potential; }

	/**
	 * @brief Returns the radius for CUSP potential in FMM (ScaFaCoS). 
	 *
	 * @return Radius for CUSP potential in FMM. 
	 */
	double isFCS_radius();

	/**
	 * @brief Returns the units for the radius of CUSP potentials in FMM (ScaFaCoS). 
	 *
	 * @return QString containing the units for the radius. 
	 */
	QString isFCS_radiusUnit();

	/**
	 * @brief Returns the limit for unrolled functions in FMM (ScaFaCoS). 
	 * 
	 * @return Limit for unrolled functions in FMM. 
	 */
	long long isFCS_unroll_limit() const { return fcs_unroll_limit; }

	/**
	 * @brief Returns the debug level for PEPC (ScaFaCoS). 
	 * 
	 * @return Debug level for PEPC. 
	 */
	int isFCS_debuglevel() const { return fcs_debuglevel; }

	/**
	 * @brief Returns the epsilon value for PEPC (ScaFaCoS). 
	 * 
	 * @return Epsilon for PEPC. 
	 */
	double isFCS_epsilon();

	/**
	 * @brief Returns the Units for the epsilon value for PEPC (ScaFaCoS). 
	 *
	 * @return QString containing the units for the epsilon value. 
	 */
	QString isFCS_epsilonUnit(); 

	/**
	 * @brief Returns the value of the load balancing switch in PEPC (ScaFaCoS). 
	 *
	 * @return Load balancing switch for PEPC. 
	 */
	int isFCS_load_balancing() const { return fcs_load_balancing; }

	/**
	 * @brief Returns NPM for PEPC (ScaFaCoS). 
	 * 
	 * @return NPM for PEPC. 
	 */
	double isFCS_npm(); 

	/**
	 * @brief Returns the Units for the NPM in PEPC (ScaFaCoS). 
	 *
	 * @return QString containing the units for the NPM. 
	 */ 
	QString isFCS_npmUnit();

	/**
	 * @brief Returns the number of walk threads per MPI rank for PEPC (ScaFaCoS). 
	 *
	 * @return Number of walk threads for PEPC. 
	 */ 
	int isFCS_num_walk_threads() const { return fcs_num_walk_threads; }

	/**
	 * @brief Returns the theta value for PEPC (ScaFaCoS). 
	 * 
	 * @return Theta value for PEPC. 
	 */
	double isFCS_theta();

	/**
	 * @brief Returns the units for the theta value in PEPC (ScaFaCoS). 
	 * 
	 * @return QString containing the units for the theta value. 
	 */
	QString isFCS_thetaUnit();

	/**
	 * @brief Returns the interpolation degree for PP3MG (ScaFaCoS).
	 *
	 * @return Interpolation degree for PP3MG. 
	 */ 
	int isFCS_degree() const { return fcs_degree; }

	/**
	 * @brief Returns the number of ghost cells for PP3MG (ScaFaCoS). 
	 *
	 * @return Number of ghost celles for PP3MG. 
	 */
	int isFCS_ghosts() const { return fcs_ghosts; }

	/**
	 * @brief Returns the number of grid cells in x-dimension for PP3MG (ScaFaCoS). 
	 *
	 * @return Number of grid cells in x-dimension for PP3MG. 
	 */ 
	int isFCS_gridsize_x() const { return fcs_gridsize_x.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the number of grid cells in y-dimension for PP3MG (ScaFaCoS). 
	 *
	 * @return Number of grid cells in y-dimension for PP3MG. 
	 */ 
	int isFCS_gridsize_y() const { return fcs_gridsize_y.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the number of grid cells in z-dimension for PP3MG (ScaFaCoS). 
	 *
	 * @return Number of grid cells in z-dimension for PP3MG. 
	 */ 
	int isFCS_gridsize_z() const { return fcs_gridsize_z.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the number of iterations for PP3MG and VMG (ScaFaCoS). 
	 *
	 * @return Number of iterations for PP3MG and VMG. 
	 */
	int isFCS_max_iterations() const { return fcs_max_iterations.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the multi-grid cycle type specifier for VMG (ScaFaCoS). 
	 * 
	 * @return Multi-grid cycle type specifier. 
	 */
	int isFCS_cycle_type() const { return fcs_cycle_type; }

	/**
	 * @brief Returns the order of the discretization scheme for PDE solver for VMG (ScaFaCoS).
	 *
	 * @return Order of the discretization scheme. 
	 */
	int isFCS_discretization_order() const { return fcs_discretization_order; }

	/**
	 * @brief Returns the interpolation order for VMG (ScaFaCoS). 
	 *
	 * @return Degree of interpolation polynomial. 
	 */
	int isFCS_interpolation_order() const { return fcs_interpolation_order.at(myLongrangeAlgo); }

	/**
	 * @brief Returns the maximum number of multi-grid levels for VMG (ScaFaCoS). 
	 *
	 * @return Maximum number of multi-grid levels. 
	 */
	int isFCS_max_level() const { return fcs_max_level; }

	/**
	 * @brief Returns the range of smearing of charges over near field cells for VMG (ScaFaCoS).
	 *
	 * @return Range of smearing of charges over near field cells for VMG. 
	 */
	int isFCS_near_field_cells() const { return fcs_near_field_cells; }

	/**
	 * @brief Returns residual threshold for solver iteration in VMG (ScaFaCoS). 
	 *
	 * @return Residual threshold for VMG. 
	 */
	double isFCS_precision();

	/**
	 * @brief Returns the unit for the residual threshold for solver iteration in VMG (ScaFaCoS). 
	 *
	 * @return QString containing the units for the residual threshold. 
	 */
	QString isFCS_precisionUnit();

	/**
	 * @brief Returns number of smoothing steps in VMG (ScaFaCoS). 
	 *
	 * @return Number of (pre- and post-) smoothing steps. 
	 */
	int isFCS_smoothing_steps() const { return fcs_smoothing_steps; }

	/**
	 * @brief Returns value of epsI for P2NFFT (ScaFaCoS). 
	 *
	 * @return Value of epsI. 
	 */
	double isFCS_epsI();

	/**
	 * @brief Returns the unit of epsI for P2NFFT (ScaFaCoS). 
	 *
	 * @return QString containing the units for the epsI value. 
	 */
	QString isFCS_epsIUnit();

	/**
	 * @brief Returns value of parameter m for P2NFFT (ScaFaCoS). 
	 *
	 * @return Value of m. 
	 */
	int isFCS_m() const { return fcs_m; }

	/**
	 * @brief Returns size of oversampled FFT grid in x-dimension for P2NFFT (ScaFaCoS). 
	 *
	 * @return Size of oversampled FFT grid in x-dimension.
	 */
	int isFCS_oversampled_gridsize_x() const {return fcs_oversampled_gridsize_x; }

	/**
	 * @brief Returns size of oversampled FFT grid in y-dimension for P2NFFT (ScaFaCoS). 
	 *
	 * @return Size of oversampled FFT grid in y-dimension.
	 */
	int isFCS_oversampled_gridsize_y() const {return fcs_oversampled_gridsize_y; }

	/**
	 * @brief Returns size of oversampled FFT grid in z-dimension for P2NFFT (ScaFaCoS). 
	 *
	 * @return Size of oversampled FFT grid in z-dimension.
	 */
	int isFCS_oversampled_gridsize_z() const {return fcs_oversampled_gridsize_z; }

	/**
	 * @brief Returns value of the parameter p used in P2NFFT (ScaFaCoS). 
	 *
	 * @return Value of parameter p.
	 */
	int isFCS_p() const {return fcs_p;}

	/**
	 * @brief Indicates the use of one number of processes, instead of cutting the domain.
	 *
	 * @return True, if Barnes-Hut or FMM are selected. 
	 */
	bool isSingleProcs() { return useSingleProcs; }

	/**
	 * @brief Indicates the number of cuts of the domain in X-direction for parallelization. 
	 *
	 * @return The number of cuts in X-direction.
	 */
	int isCutsX() { return numberOfProcs[0]; }

	/**
	 * @brief Indicates the number of cuts of the domain in Y-direction for parallelization. 
	 *
	 * @return The number of cuts in Y-direction.
	 */
	int isCutsY() { return numberOfProcs[1]; }
	
	/**
	 * @brief Indicates the number of cuts of the domain in Y-direction for parallelization. 
	 *
	 * @return The number of cuts in Y-direction.
	 */
	int isCutsZ() { return numberOfProcs[2]; }

	/**
	 * @brief Returns the number of processes used in the case of Barnes-Hut or FMM. 
	 *
	 * @return Number of processes.
	 */
	int isProcsAtAll() { return procsAtAll; }

	/**
	 * @brief Returns the total number of processes used (except for Barnes-Hut or FMM).
	 *
	 * @return Total number of processes. 
	 */
	long isProcessorUsage() const;
	
	/**
	 * @brief Returns the list of available loadbalancing functions. 
	 *
	 * @return A List containing all available loadbalancing functions. 
	 */
	QStringList* isLoadBalancingFunctionList() { return loadBalancingFunctionList; }

	/**
	 * @brief Returns the name of the actual selected loadbalancing function. 
	 *
	 * @return A string containing the name of the current loadbalancing function. 
	 */
	QString isLoadBalancingFunction() { return *loadBalancingFunction; }

	/**
	 * @brief Returns the starttime for loadbalancing. 
	 *
	 * @return Starttime for loadbalancing. 
	 */
	double isLoadbal_t();

	/**
	 * @brief Returns the starttime's unit for loadbalancing. 
	 *
	 * @return A string indicating the unit for the loadbalancing starttime. 
	 */
	QString isLoadbal_tUnit();

	/**
	 * @brief Returns the timestep width for loadbalancing. 
	 * 
	 * @return Timedelta for loadbalancing. 
	 */ 
	double isLoadbal_delta();

	/**
	 * @brief Returns the unit of the timedelta value for loadbalancing. 
	 *
	 * @return A string indicating the unit for LB timedelta. 
	 */
	QString isLoadbal_deltaUnit();
	
	/**
	 * @brief Returns the steps for loadbalancing. 
	 *
	 * @return The steps for loadbalancing. 
	 */ 
	int isLoadbal_step() { return loadbal_step; }

public slots:
	
	/**
	 * @brief Resets all parameter values to the defaults (this is not reversible).
	 *        Note that the parameters of all longrange algorithms are reset!
	 */
	void clear();

	/**
	 * @brief Selects the longrange solvers of this SolParallel_Data object and enables the corresponding parameters. 
	 * 
	 * @param algo Index of the longrange algorithm corresponding to #longrange_algorithms.
	 */
	void setLongrangeAlgo( int algo );

	/**
	 * @brief Sets the linked cell r_cut parameter and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the linked cell r_cut.
	 */
	void setLinkedCellR_Cut( double valueIn );

	/**
	 * @brief Sets the coulomb force constant and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the coulomb force constant. 
	 */
	void setEpsilon0inv( double valueIn );

	/**
	 * @brief Sets the r_cut of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the r_cut parameter of the current longrange solver. 
	 */
	void setR_Cut( double valueIn );

	/**
	 * @brief Sets the r_l value of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the r_l parameter of the current longrange solver. 
	 */	 
	void setR_L( double valueIn );

	/**
	 * @brief Sets the splitting coefficient of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the splitting coefficient of the current longrange solver. 
	 */
	void setSplittingCoefficientG( double valueIn );

	/**
	 * @brief Sets the MAP value of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the MAP parameter of the current longrange solver. 
	 */
	void setMAP( double valueIn );

	/**
	 * @brief Sets the celleratio of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the cellratio of the current longrange solver. 
	 */
	void setCellratio( int valueIn );
	
	/**
	 * @brief Sets the interpolation degree of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the interpolation degree of the current longrange solver. 
	 */
	void setInterpolationDegree( int valueIn );

	/**
	 * @brief Sets the maximum tree level of the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the maximum tree level of the current longrange solver. 
	 */
	void setMaxTreeLevel( int valueIn );

	/**
	 * @brief Sets the poisson solver for the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn Name of the poisson solver for the current longrange solver. 
	 */
	void setPoisson_solver( const QString &valueIn );

	/**
	 * @brief Sets the tolerance (for ScaFaCoS solver only) for the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the tolerance. 
	 */
	void setFCS_tolerance(double valueIn);

	/**
	 * @brief Sets the tolerance type (for ScaFaCoS solver only) for the current longrange algorithm and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the tolerance type (see ScaFaCoS manual). 
	 */
	void setFCS_tolerance_type(int valueIn);

	/**
	 * @brief Sets the splitting coefficient for Ewald method/P2NFFT (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the splitting coefficient. 
	 */
	void setFCS_splittingCoefficientAlpha(double valueIn);

	/**
	 * @brief Sets the K-space cut-off for the Ewald method (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the K-space cut-off. 
	 */
	void setFCS_kmax(int valueIn);

	/**
	 * @brief Sets the maximum K-space cut-off for the Ewald method (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the maximum K-space cut-off. 
	 */
	void setFCS_maxkmax(int valueIn); 

	/**
	 * @brief Sets the number of periodic images used in x-dimension for the direct method  (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the number of periodic images in x-dimension.
	 */
	void setFCS_periodic_images_x(int valueIn);

	/**
	 * @brief Sets the number of periodic images used in y-dimension for the direct method  (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the number of periodic images in y-dimension.
	 */
	void setFCS_periodic_images_y(int valueIn);

	/**
	 * @brief Sets the number of periodic images used in z-dimension for the direct method  (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the number of periodic images in z-dimension.
	 */
	void setFCS_periodic_images_z(int valueIn);
	
	/**
	 * @brief Sets the load balancing routine for FMM (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn Index for the load balancing routine. 
	 */
	void setFCS_balanceload(int valueIn);

	/**
	 * @brief Sets the dipole correction for FMM and PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the dipole correction. 
	 */
	void setFCS_dipole_correction(int valueIn);

	/**
	 * @brief Sets the maximum tree depth for FMM (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the maximum tree depth. 
	 */
	void setFCS_maxdepth(long long valueIn);

	/**
	 * @brief Sets the potential used for FMM (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New index for the potential. 
	 */
	void setFCS_potential(int valueIn);

	/**
	 * @brief Sets the radius for CUSP potential in FMM (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the radius.
	 */
	void setFCS_radius(double valueIn);

	/**
	 * @brief Sets the limit for unrolled functions in FMM (ScaFaCoS) and emits a hasChanged() signal.
	 *
	 * @param valueIn New value for the limit for unrolled functions. 
	 */
	void setFCS_unroll_limit(long long valueIn);

	/**
	 * @brief Sets the debug level for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the debug level. 
	 */
	void setFCS_debuglevel(int valueIn);

	/**
	 * @brief Sets the epsilon value for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New epsilon value. 
	 */
	void setFCS_epsilon(double  valueIn);

	/**
	 * @brief Sets the load balancing switch for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the load balancing switch.
	 */
	void setFCS_load_balancing(int valueIn);

	/**
	 * @brief Sets the NPM value for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New NPM value.
	 */
	void setFCS_npm(double valueIn);

	/**
	 * @brief Sets the number of walk threads per MPI rank for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of walk threads. 
	 */
	void setFCS_num_walk_threads(int valueIn);

	/**
	 * @brief Sets the theta value for PEPC (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New theta value. 
	 */
	void setFCS_theta(double valueIn);

	/**
	 * @brief Sets the interpolation degree for PP3MG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the interpolation degree. 
	 */
	void setFCS_degree(int valueIn);

	/**
	 * @brief Sets the number of ghost cells for PP3MG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of ghost cells. 
	 */
	void setFCS_ghosts(int valueIn);

	/**
	 * @brief Sets the number of grid cells in x-dimension for PP3MG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of grid cells in x-dimension. 
	 */
	void setFCS_gridsize_x(int valueIn);

	/**
	 * @brief Sets the number of grid cells in y-dimension for PP3MG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of grid cells in y-dimension. 
	 */
	void setFCS_gridsize_y(int valueIn);

	/**
	 * @brief Sets the number of grid cells in z-dimension for PP3MG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of grid cells in z-dimension. 
	 */
	void setFCS_gridsize_z(int valueIn);

	/**
	 * @brief Sets the maximum number of iterations for PP3MG and VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New maximum number of iterations. 
	 */
	void setFCS_max_iterations(int valueIn);

	/**
	 * @brief Sets the multi-grid cycle type specifier for VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the multi-grid cycle type specifier. 
	 */
	void setFCS_cycle_type(int valueIn);

	/**
	 * @brief Sets the discretization order scheme in VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the discretization order scheme. 
	 */
	void setFCS_discretization_order(int valueIn);

	/**
	 * @brief Sets the degree of interpolation polynomial in VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the interpolation degree. 
	 */
	void setFCS_interpolation_order(int valueIn);

	/**
	 * @brief Sets the maximum number of multi-grid levels for VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the maximum number of multi-grid levels. 
	 */
	void setFCS_max_level(int valueIn);

	/**
	 * @brief Sets the range of smearing of charges over near field cells in VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the range of smearing of charges over near field cells. 
	 */
	void setFCS_near_field_cells(int valueIn);

	/**
	 * @brief Sets the threshold for residual for solver iteration in VMG (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for the threshold. 
	 */
	void setFCS_precision(double valueIn);

	/**
	 * @brief Sets the number of smoothing steps in VMG (ScaFaCoS) and emits a hasChanged() signal.
	 *
	 * @param valueIn New value for the number of smoothing steps. 	 
	 */
	void setFCS_smoothing_steps(int valueIn);

	/**
	 * @brief Sets the value of epsI used in P2NFFT (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for epsI parameter. 
	 */
	void setFCS_epsI(double valueIn);

	/**
	 * @brief Sets the value of parameter m used in P2NFFT (ScaFaCoS) and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value for parameter m. 
	 */
	void setFCS_m(int valueIn);

	/**
	 * @brief Sets the size of the oversampled FFT grid in x-dimension for P2NFFT (ScaFaCoS) and emits  a hasChanged() signal. 
	 *
	 * @param valueIn New value. 
	 */
	void setFCS_oversampled_gridsize_x(int valueIn);

	/**
	 * @brief Sets the size of the oversampled FFT grid in y-dimension for P2NFFT (ScaFaCoS) and emits  a hasChanged() signal. 
	 *
	 * @param valueIn New value. 
	 */
	void setFCS_oversampled_gridsize_y(int valueIn);

	/**
	 * @brief Sets the size of the oversampled FFT grid in z-dimension for P2NFFT (ScaFaCoS) and emits  a hasChanged() signal. 
	 *
	 * @param valueIn New value. 
	 */
	void setFCS_oversampled_gridsize_z(int valueIn);

	/**
	 * @brief Sets the value of the parameter p used in P2NFFT (ScaFaCoS) and emits  a hasChanged() signal. 
	 *
	 * @param valueIn New value for the parameter p. 
	 */
	void setFCS_p(int valueIn);

	/**
	 * @brief Sets the number of processes in X-direction and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of processes in X-direction. 
	 */
	void setCutsX( int cutsIn);

	/**
	 * @brief Sets the number of processes in Y-direction and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of processes in Y-direction. 
	 */
	void setCutsY( int cutsIn);

	/**
	 * @brief Sets the number of processes in Z-direction and emits a hasChanged() signal. 
	 *
	 * @param valueIn New number of processes in Z-direction. 
	 */
	void setCutsZ( int cutsIn);

	/**
	 * @brief Sets the number of processes in the case of Barnes-Hut or FMM and emits a hasChanged() signal. 
	 *
	 * @param valueIn Number of processes. 
	 */
	void setProcsAtAll( int procsIn );

	/**
	 * @brief Sets the load balancing function and emits a hasChanged() signal. 
	 *
	 * @param valueIn String containing the name of the load balancing function. 
	 */
	void setLoadBalancingFunction( const QString &valueIn );

	/**
	 * @brief Sets the starttime for load balancing and emits a hasChanged() signal. 
	 *
	 * @param valueIn New starttime for load balancing. 
	 */
	void setLoadbal_t( double valueIn );

	/**
	 * @brief Sets the timedelta for load balancing and emits a hasChanged() signal. 
	 *
	 * @param valueIn New value of timedelta for load balancing. 
	 */
	void setLoadbal_delta( double valueIn );

	/**
	 * @brief Sets the steps for load balancing and emits a hasChanged() signal. 
	 *
	 * @param valueIn Steps for load balancing. 
	 */
	void setLoadbal_step( int valueIn );
     
signals:
	
	/*
	 * Signals used to enable parameters for longrange solvers and parallelization. 
	 */
	void enablepoisson_solver(bool enable);
	void enableR_Cut(bool enable);
	void enableR_L(bool enable);
	void enableSplittingCoefficientG(bool enable);
	void enableMAP(bool enable);
	void enableMaxTreeLevel(bool enable);
	void enableInterpolationDegree(bool enable);
	void enableCellratio(bool enable);
	void enableFCS_tolerance(bool enable);
	void enableFCS_tolerance_type(bool enable);
	void enableFCS_periodic_images_x(bool enable);
	void enableFCS_periodic_images_y(bool enable);
	void enableFCS_periodic_images_z(bool enable);
	void enableFCS_splittingCoefficientAlpha(bool enable);
	void enableFCS_kmax(bool enable);
	void enableFCS_maxkmax(bool enable);
	void enableFCS_balanceload(bool enable);
	void enableFCS_dipole_correction(bool enable);
	void enableFCS_maxdepth(bool enable);
	void enableFCS_potential(bool enable);
	void enableFCS_radius(bool enable);
	void enableFCS_unroll_limit(bool enable);
	void enableFCS_degree(bool enable);
	void enableFCS_ghosts(bool enable);
	void enableFCS_gridsize_x(bool enable);
	void enableFCS_gridsize_y(bool enable);
	void enableFCS_gridsize_z(bool enable);
	void enableFCS_max_iterations(bool enable);
	void enableFCS_debuglevel(bool enable);
	void enableFCS_epsilon(bool enable);
	void enableFCS_load_balancing(bool enable);
	void enableFCS_npm(bool enable);
	void enableFCS_num_walk_threads(bool enable);
	void enableFCS_theta(bool enable);
	void enableFCS_cycle_type(bool enable);
	void enableFCS_discretization_order(bool enable);
	void enableFCS_interpolation_order(bool enable);
	void enableFCS_max_level(bool enable);
	void enableFCS_near_field_cells(bool enable);
	void enableFCS_precision(bool enable);
	void enableFCS_smoothing_steps(bool enable);
	void enableFCS_epsI(bool enable);
	void enableFCS_m(bool enable);
	void enableFCS_p(bool enable);
	void enableFCS_oversampled_gridsize_x(bool enable);
	void enableFCS_oversampled_gridsize_y(bool enable);
	void enableFCS_oversampled_gridsize_z(bool enable);
	void enableLoadBalValues_dyn(bool enable);
	void enableLoadBalValues_opt(bool enable);

     
     
protected:
     
protected slots:
	/*
	 * @brief Emits a hasChanged() signal if the given parameter equals cAll. 
	 *
	 * @param whatShouldISet A whoIsWho metadatatype. 
	 */
	void setAnyOption( whoIsWho whatShouldISet );

private:
	 
	/// Parent simulationParameterData class. 
	simulationParameterData* myParent;
	
	/// Gui configuration for this SolParallel_Data object. 
	TremoloGUIConfig* myConfig;
	
	/// List of keywords for the 'Solver and Parallelization' section of the .parameters file.
	QStringList* myKeywordList;
	
	/// Linked Cell r_cut.
	StandardisedDouble linkedCellR_Cut;

	/// Coulomb force constant.
	StandardisedDouble epsilon0inv;

	/// The current longrange solver of this SolParallel_Data class. 
	longrange_algorithms myLongrangeAlgo;
	 
	/**
	 * @brief The parameters for the longrange solvers.
	 * 
	 * @note All parameters that are used by multiple solvers are stored in vectors of length size_count. 
	 *       All parameters that are used by only a single solver are stored as separate variables. 
	 */

	/// Longrange solver's r_cut. @todo consistency check. 
	std::vector<StandardisedDouble> r_cut;
	
	/// N^2Spline r_l.
	StandardisedDouble r_l;
	 
	/// SPME splitting coefficient.
	StandardisedDouble splittingCoefficientG;

	/// FMM MAP.
	StandardisedDouble MAP;

	/// SPME cellratio. 
	int cellratio;

	/// Longrange solver's interpolation degree.
	std::vector<int> interpolationDegree;

	/// Longrange solver's maximum tree level. 
	int maxTreeLevel;

	/// Poisson solver of longrange solver. @todo consistency check. 
	QString* poissonsolver;

	/// List of available poisson solvers (currently "FFT" and "multigrid").
	QStringList* poissonsolverList;

	/// List containing the strings used in the parameter file for all available poisson solvers.
	QStringList* poissonsolverFileList;

	/**
	 * @brief Tolerance (ScaFaCoS).
	 *
	 * @todo Determine correct data type!
	 */
	std::vector<StandardisedDouble> fcs_tolerance;

	/**
	 * @brief Tolerance type (ScaFaCoS).
	 */
	std::vector<int> fcs_tolerance_type;

	/**
	 * @brief Splitting parameter /f$alpha/f$ for Ewald method and P2NFFT (ScaFaCoS). 
	 *
	 * @todo Determine correct data type! Is this the same parameter as in SPME (splittingCoefficientG)?
	 */
	std::vector<StandardisedDouble> fcs_splittingCoefficientAlpha;

	/**
	 * @brief K-space cut-off for Ewald method (ScaFaCoS).
	 *
	 * @todo Determine correct data type!
	 */
	int fcs_kmax;

	/**
	 * @brief Maximum K-space cut-off for Ewald method (ScaFaCos).
	 *
	 * @todo Determine correct data type!
	 */
	int fcs_maxkmax;

	/// Number of periodic images used in x-dimension (ScaFaCoS). 
	int fcs_periodic_images_x;

	/// Number of periodic images used in y-dimension (ScaFaCoS). 
	int fcs_periodic_images_y;

	/// Number of periodic images used in z-dimension (ScaFaCoS). 
	int fcs_periodic_images_z;

	/// Load balancing routine for FMM (ScaFaCoS).  
	int fcs_balanceload;

	/// Dipole correction for FMM and PEPC (ScaFaCoS). 
	std::vector<int> fcs_dipole_correction;

	/// Maximum tree depth for FMM (ScaFaCoS). 
	long long fcs_maxdepth;

	/// Potential used for FMM (ScaFaCoS)
	int fcs_potential;

	/// Radius for CUSP potential in FMM (ScaFaCoS). 
	StandardisedDouble fcs_radius;

	/// Limit for unrolled functions in FMM (ScaFaCoS).
	long long fcs_unroll_limit;

	/// Debug level for PEPC (ScaFaCoS). 
	int fcs_debuglevel;

	/// Epsilon for PEPC (ScaFaCoS). 
	StandardisedDouble fcs_epsilon;

	/// Load balancing switch for PEPC (ScaFaCoS). 
	int fcs_load_balancing;

	/// NPM for PEPC (ScaFaCoS).
	StandardisedDouble fcs_npm;

	/// Number of walk threads per MPI rank for PEPC (ScaFaCoS). 
	int fcs_num_walk_threads;

	/// Theta for PEPC (ScaFaCoS). 
	StandardisedDouble fcs_theta;

	/// Interpolation degree for PP3MG (ScaFaCoS). 
	int fcs_degree;

	/// Number of ghost cells for PP3MG (ScaFaCoS). 
	int fcs_ghosts;

	/// Number of grid cells in x-dimension for PP3MG (ScaFaCoS). 
	std::vector<int> fcs_gridsize_x;

	/// Number of grid cells in y-dimension for PP3MG (ScaFaCoS). 
	std::vector<int> fcs_gridsize_y;

	/// Number of grid cells in z-dimension for PP3MG (ScaFaCoS). 
	std::vector<int> fcs_gridsize_z;

	/// Maximum number of iterations for PP3MG and VMG (ScaFaCoS). 
	std::vector<int> fcs_max_iterations;

	/// Multi-grid cycle type specifier for VMG (ScaFaCoS). 
	int fcs_cycle_type;

	/// Order of discretization scheme for PDE solver for VMG (ScaFaCoS). 
	int fcs_discretization_order;

	/// Degree of interpolation polynomial for VMG (ScaFaCoS). 
	std::vector<int> fcs_interpolation_order;

	/// Maximum number of multi-grid levels for VMG (ScaFaCoS). 
	int fcs_max_level;

	/// Range of smearing of charges over near field cells for VMG (ScaFaCoS). 
	int fcs_near_field_cells;

	/// Threshold for residual for solver iteration for VMG (ScaFaCoS). 
	StandardisedDouble fcs_precision;

	/// Number of pre- and post-smoothing steps for VMG (ScaFaCoS). 
	int fcs_smoothing_steps; 

	/// epsI parameter for P2NFFT (ScaFaCoS). 
	StandardisedDouble fcs_epsI;

	/// Parameter m for P2NFFT (ScaFaCoS).
	int fcs_m;

	/// Size of oversampled FFT grid in x-dimension used in P2NFFT (ScaFaCoS). 
	int fcs_oversampled_gridsize_x;

	/// Size of oversampled FFT grid in y-dimension used in P2NFFT (ScaFaCoS).  
	int fcs_oversampled_gridsize_y;

	/// Size of oversampled FFT grid in z-dimension used in P2NFFT (ScaFaCoS). 
	int fcs_oversampled_gridsize_z;

	/// Parameter p used in P2NFFT (ScaFaCoS.)
	int fcs_p;

	/*
	 * All parameters related with the parallelization settings. 
	 */
	/// Load balancing function for parallelization
	QString* loadBalancingFunction;

	/// List of all available load balancing functions (currently "OFF", "LINEAR" and "QUADRATIC"). 
	QStringList* loadBalancingFunctionList;

	/// Number of processes in every coordinate direction. 
	int numberOfProcs[3];

	/// Overall number of used processes.
	int procsAtAll;

	/// Indicates if single processes are used (true except for Barnes-Hut and FMM).
	bool useSingleProcs;

	/// Starttime for load balancing. 
	StandardisedDouble loadbal_t;

	/// Timedelta for load balancing.
	StandardisedDouble loadbal_delta;

	/// Steps for load balancing. 
	int loadbal_step;

private slots:
	/**
	 * @brief Sets the type of integration (optimization or dynamics)
	 */
	void change_Integration( bool dynamo );
};
 
#endif
