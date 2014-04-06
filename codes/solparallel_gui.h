/***************************************************************************
                          solparallel_gui.h  -  description
                             -------------------
    begin                : Fri Jun 04 2004
			   Written by Daniel Bloemer, Tim Golla
    email                : bloemer@iam.uni-bonn.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   Tremolo-X Project                                                     *
 *   All rights reserved                                                   *
 *   For details please see the LICENSE file                               *
 *                                                                         *
 ***************************************************************************/

#ifndef HAVE_SOLPARALLEL_GUI
#define HAVE_SOLPARALLEL_GUI

// QT-includes
#include <QWidget>
#include <QGridLayout>
#include <Q3VBoxLayout>
#include <QLabel>
#include <Q3Frame>

class SolParallel_Data;
class TremoloGUIConfig;
class Q3ButtonGroup;
class QComboBoxWR;
class Q3Frame;
class QGridLayout;
class QLabel;
class QLCDNumber;
class QLineEdit;
class QRadioButton;
class QSpinBox;
class Q3WidgetStack;
class Q3VBoxLayout;


/**
 * @brief This class provides the GUI for solver and parallelization
 * @author Daniel Bloemer <bloemer@iam.uni-bonn.de>
 * @author Tim Golla <golla@ins.uni-bonn.de>
 *
 * @warning Some of the used Qt classes may be outdated. 
 */
class SolParallel_GUI : public QWidget
{
	Q_OBJECT

public:
	 
	/// Constructor.
	SolParallel_GUI( SolParallel_Data* DataPart, TremoloGUIConfig* configIn, QWidget* parent = 0, const char* name = 0, Qt::WFlags f = 0 );

	/// Destructor.
	~SolParallel_GUI();

public slots:
     
signals:
     
protected:
     
protected slots:

	/**
	 * @brief Slot responsible for language changes. 
	 */
	virtual void languageChange();

    /**
	 * @brief Slot updating all data related to the SolParallel_Data class. 
	 */
	void updateGUI_data();

	/**
	 * @brief Slot updating all data related to the SolParallel_Data class and additionally the units of all parameters. 
	 */ 
	void updateGUI_all();

	void lCR_CLErp();
	void epsilon0inv_rp();
	void R_CLErp();
	void R_LLErp();
	void sCGLErp();
	void mapLErp();
	void cellratioLErp();
	void interpolationLErp();
	void maxTreeLErp();
	void fcstoleranceLErp();
	void fcstoleranceTypeLErp();
	void fcssCAlphaLErp();
	void fcspImagesXLErp();
	void fcspImagesYLErp();
	void fcspImagesZLErp();
	void fcskmaxLErp();
	void fcsmaxkmaxLErp();
	void fcsbalanceloadLErp();
	void fcsdipoleCLErp();
	void fcsmaxdepthLErp();
	void fcspotentialLErp();
	void fcsradiusLErp();
	void fcsuLimitLErp();
	void fcsdegreeLErp();
	void fcsghostsLErp();
	void fcsgSizeXLErp();
	void fcsgSizeYLErp();
	void fcsgSizeZLErp();
	void fcsmIterationsLErp();
	void fcsdebuglevelLErp();
	void fcsepsilonLErp();
	void fcsloadbalancingLErp();
	void fcsnpmLErp();
	void fcsnwthreadsLErp();
	void fcsthetaLErp();
	void fcscycletypeLErp();
	void fcsdisOrderLErp();
	void fcsintOrderLErp();
	void fcsmaxLevelLErp();
	void fcsnfcellsLErp();
	void fcsprecisionLErp();
	void fcssStepsLErp();
	void fcsepsILErp();
	void fcsmLErp();
	void fcspLErp();
	void fcsoversampledGsXLErp();
	void fcsoversampledGsYLErp();
	void fcsoversampledGsZLErp();
	void lb_tLErp();
	void lb_deltaLErp();
	void lb_stepLErp();

private:
	
	/// Central layout.
	Q3VBoxLayout* centralLayout;

	/// Buttongroup for the longrange solvers. 
	Q3ButtonGroup* LongrangeButtonGroup;

	/// Layout for the longrange solver's buttongroup. 
	QGridLayout* LongrangeButtonGroupLayout;

	// Radiobuttons for all longrage solvers. 
	QRadioButton* OFFRadioButton;
	QRadioButton* N_2RadioButton;
	QRadioButton* N_2SplineRadioButton;
	QRadioButton* SPMERadioButton;
	QRadioButton* FMMRadioButton;
	QRadioButton* FCS_DirectRadioButton;
	QRadioButton* FCS_EwaldRadioButton;
	QRadioButton* FCS_FMMRadioButton;
	QRadioButton* FCS_PEPCRadioButton;
	QRadioButton* FCS_PP3MGRadioButton;
	QRadioButton* FCS_VMGRadioButton;
	QRadioButton* FCS_P2NFFTRadioButton;

	/// Buttongroup for parameter widgets. 
	Q3ButtonGroup* DataButtonGroup;

	/// Layout for the parameter's buttongroup. 
	QGridLayout* DataButtonGroupLayout;

	/// Label for linked cell r_cut. 
	QLabel* linkedCellR_CutLabel;

	/// LineEdit for linked cell r_cut. 
	QLineEdit* linkedCellR_CutLineEdit;

	/// Label for coulomb force constant. 
	QLabel* epsilon0invLabel;

	/// LineEdit for coulomb force constant. 
	QLineEdit* epsilon0invLineEdit;

	/// Label for r_cut parameter. 
	QLabel* r_curLabel;

	/// LineEdit for r_cut parameter. 
	QLineEdit* r_cutLineEdit;

	/// Label for r_l parameter. 
	QLabel* r_lLabel;

	/// LineEdit for r_l parameter. 
	QLineEdit* r_lLineEdit;

	/// Label for splitting coefficient. 
	QLabel* splittingCoeffizientGLabel;

	/// LineEdit for splitting coefficient. 
	QLineEdit* splittingCoeffizientGLineEdit;

	/// Label for MAP parameter. 
	QLabel* MapLabel;

	/// LineEdit for MAP parameter. 
	QLineEdit* MapLineEdit;

	/// Label for cellratio. 
	QLabel* cellratioLabel;

	/// LineEdit for cellratio. 
	QLineEdit* cellratioLineEdit;

	/// Label for interpolation degree. 
	QLabel* interpolationLabel;

	/// LineEdit for interpolation degree. 
	QLineEdit* interpolationLineEdit;

	/// Label for maximum tree level. 
	QLabel* maxTreeLabel;

	/// LineEdit for maximum tree level. 
	QLineEdit* maxTreeLineEdit;

	/// Label for poisson solver. 
	QLabel* poisson_solverLabel;

	/// ComboBox for poisson solver.
	QComboBoxWR* poisson_solverComboBox;

	/// Label for splitting coefficient alpha (ScaFaCoS). 
	QLabel* fcsSplittingCoefficientAlphaLabel;

	/// LineEdit for splitting coefficient alpha (ScaFaCoS). 
	QLineEdit* fcsSplittingCoefficientAlphaLineEdit;

	/// Label for tolerance (ScaFaCoS). 
	QLabel* fcsToleranceLabel;

	/// LineEdit for tolerance (ScaFaCoS). 
	QLineEdit* fcsToleranceLineEdit;

	/// Label for tolerance_type (ScaFaCoS). 
	QLabel* fcsToleranceTypeLabel;
	
	/// LineEdit for tolerance_type (ScaFaCoS). 
	QLineEdit* fcsToleranceTypeLineEdit;

	/// Label for periodic_images (ScaFaCoS).
	QLabel* fcsPeriodicImagesLabel;
	
	/// LineEdit for periodic_images_x (ScaFaCoS). 
	QLineEdit* fcsPeriodicImagesXLineEdit;

	/// LineEdit for periodic_images_y (ScaFaCoS). 
	QLineEdit* fcsPeriodicImagesYLineEdit;

	/// LineEdit for periodic_images_z (ScaFaCoS). 
	QLineEdit* fcsPeriodicImagesZLineEdit;

	/// Label for kmax (ScaFaCoS). 
	QLabel* fcsKmaxLabel;
	
	/// LineEdit for kmax (ScaFaCoS). 
	QLineEdit* fcsKmaxLineEdit;

	/// Label for maxkmax (ScaFaCoS). 
	QLabel* fcsMaxkmaxLabel;

	/// LineEdit for maxkmax (ScaFaCoS).
	QLineEdit* fcsMaxkmaxLineEdit;

	/// Label for balanceload (ScaFaCoS). 
	QLabel* fcsBalanceloadLabel;

	/// LineEdit for balanceload (ScaFaCoS). 
	QLineEdit* fcsBalanceloadLineEdit;

	/// Label for dipole_correction (ScaFaCoS). 
	QLabel* fcsDipoleCorrectionLabel;
	
	/// LineEdit for dipole_correction (ScaFaCoS). 
	QLineEdit* fcsDipoleCorrectionLineEdit;

	/// Label for maxdepth (ScaFaCoS). 
	QLabel* fcsMaxdepthLabel;

	/// LineEdit for maxdepth (ScaFaCoS). 
	QLineEdit* fcsMaxdepthLineEdit;

	/// Label for potential (ScaFaCoS). 
	QLabel* fcsPotentialLabel;

	/// LineEdit for potential (ScaFaCoS). 
	QLineEdit* fcsPotentialLineEdit;

	/// Label for radius (ScaFaCoS). 
	QLabel* fcsRadiusLabel;
	
	/// LineEdit for radius (ScaFaCoS). 
	QLineEdit* fcsRadiusLineEdit;

	/// Label for unroll_limit (ScaFaCoS). 
	QLabel* fcsUnrollLimitLabel;

	/// LineEdit for unroll_limit (ScaFaCoS). 
	QLineEdit* fcsUnrollLimitLineEdit;

	/// Label for degree (ScaFaCoS). 
	QLabel* fcsDegreeLabel;

	/// LineEdit for degree (ScaFaCoS). 
	QLineEdit* fcsDegreeLineEdit;

	/// Label for ghosts (ScaFaCoS). 
	QLabel* fcsGhostsLabel;

	/// LineEdit for ghosts (ScaFaCoS). 
	QLineEdit* fcsGhostsLineEdit;

	/// Label for gridsize (ScaFaCoS). 
	QLabel* fcsGridsizeLabel;

	/// LineEdit for gridsize_x (ScaFaCoS). 
	QLineEdit* fcsGridsizeXLineEdit;

	/// LineEdit for gridsize_y (ScaFaCoS). 
	QLineEdit* fcsGridsizeYLineEdit;

	/// LineEdit for gridsize_z (ScaFaCoS). 
	QLineEdit* fcsGridsizeZLineEdit;

	/// Label for max_iterations (ScaFaCoS). 
	QLabel* fcsMaxIterationsLabel;

	/// LineEdit for max_iterations (ScaFaCoS). 
	QLineEdit* fcsMaxIterationsLineEdit;

	/// Label for debuglevel (ScaFaCoS). 
	QLabel* fcsDebuglevelLabel;

	/// LineEdit for debuglevel (ScaFaCoS). 
	QLineEdit* fcsDebuglevelLineEdit;

	/// Label for epsilon (ScaFaCoS). 
	QLabel* fcsEpsilonLabel;

	/// LineEdit for epsilon (ScaFaCoS). 
	QLineEdit* fcsEpsilonLineEdit;

	/// Label for load_balancing (ScaFaCoS). 
	QLabel* fcsLoadBalancingLabel;

	/// LineEdit for load_balancing (ScaFaCoS). 
	QLineEdit* fcsLoadBalancingLineEdit;

	/// Label for npm (ScaFaCoS). 
	QLabel* fcsNpmLabel;

	/// LineEdit for npm (ScaFaCoS). 
	QLineEdit* fcsNpmLineEdit;

	/// Label for num_walk_threads (ScaFaCoS). 
	QLabel* fcsNumWalkThreadsLabel;

	/// LineEdit for num_walk_threads (ScaFaCoS). 
	QLineEdit* fcsNumWalkThreadsLineEdit;

	/// Label for theta (ScaFaCoS). 
	QLabel* fcsThetaLabel;

	/// LineEdit for theta (ScaFaCoS). 
	QLineEdit* fcsThetaLineEdit;

	/// Label for cycle_type (ScaFaCoS). 
	QLabel* fcsCycleTypeLabel;

	/// LineEdit for cycle_type (ScaFaCoS). 
	QLineEdit* fcsCycleTypeLineEdit;

	/// Label for discretization_order (ScaFaCoS). 
	QLabel* fcsDiscretizationOrderLabel;

	/// LineEdit for discretization_order (ScaFaCoS). 
	QLineEdit* fcsDiscretizationOrderLineEdit;

	/// Label for interpolation_order (ScaFaCoS). 
	QLabel* fcsInterpolationOrderLabel;

	/// LineEdit for interpolation_order (ScaFaCoS). 
	QLineEdit* fcsInterpolationOrderLineEdit;

	/// Label for max_level (ScaFaCoS). 
	QLabel* fcsMaxLevelLabel;

	/// LineEdit for max_level (ScaFaCoS). 
	QLineEdit* fcsMaxLevelLineEdit;

	/// Label for near_field_cells (ScaFaCoS). 
	QLabel* fcsNearFieldCellsLabel;

	/// LineEdit for near_field_cells (ScaFaCoS). 
	QLineEdit* fcsNearFieldCellsLineEdit;

	/// Label for precision (ScaFaCoS). 
	QLabel* fcsPrecisionLabel;
	
	/// LineEdit for precision (ScaFaCoS). 
	QLineEdit* fcsPrecisionLineEdit;

	/// Label for smoothing_steps (ScaFaCoS). 
	QLabel* fcsSmoothingStepsLabel;

	/// LineEdit for smoothing_steps (ScaFaCoS). 
	QLineEdit* fcsSmoothingStepsLineEdit;

	/// Label for epsI (ScaFaCoS). 
	QLabel* fcsEpsILabel;

	/// LineEdit for epsI (ScaFaCoS). 
	QLineEdit* fcsEpsILineEdit;

	/// Label for m (ScaFaCoS). 
	QLabel* fcsMLabel;

	/// LineEdit for m (ScaFaCoS). 
	QLineEdit* fcsMLineEdit;

	/// Label for p (ScaFaCoS). 
	QLabel* fcsPLabel;

	/// LineEdit for p (ScaFaCoS). 
	QLineEdit* fcsPLineEdit;

	/// Label for oversampled_gridsize (ScaFaCoS). 
	QLabel* fcsOversampledGridsizeLabel;

	/// LineEdit for oversampled_gridsize_x (ScaFaCoS). 
	QLineEdit* fcsOversampledGridsizeXLineEdit;

	/// LineEdit for oversampled_gridsize_y (ScaFaCoS). 
	QLineEdit* fcsOversampledGridsizeYLineEdit;

	/// LineEdit for oversampled_gridsize_z (ScaFaCoS). 
	QLineEdit* fcsOversampledGridsizeZLineEdit;

	/// Buttongroup for parallelization widgets. 
	Q3ButtonGroup* parallelizerButtonGroup;

	/// Layout for the parallelization buttongroup. 
	QGridLayout* parallelizerButtonGroupLayout;

	/// Label for the load balancing function. 
	QLabel* loadBalancingfunctionLabel;
	
	/// Label for cuts in X-direction. 
	QLabel* cutsXLabel;

	/// Label for cuts in Y-direction. 
	QLabel* cutsYLabel;

	/// Label for cuts in Z-direction.
	QLabel* cutsZLabel;

	/// Label for used processes.
	QLabel* usedProcsLabel;

	/// SpinBox for cuts in X-direction.
	QSpinBox* cutsXSpinBox;

	/// SpinBox for cuts in Y-direction.
	QSpinBox* cutsYSpinBox;

	/// SpinBox for cuts in Z-direction. 
	QSpinBox* cutsZSpinBox;

	/// WidgetStack used to hide and raise the "Used processes" SpinBox and the "Used processes" value label.
	Q3WidgetStack* usedProcsWidgetStack;

	/// SpinBox for used processes. 
	QSpinBox* usedProcsSpinBox;

	/// Label indicating the number of used processes. 
	QLabel* usedProcsValueLabel;

	/// ComboBox for the load balancing function. 
	QComboBoxWR* loadBalancingfunctionComboBox;
	
	/// Label for load balancing starttime. 
	QLabel* loadBal_tLabel;

	/// Label for load balancing timedelta.
	QLabel* loadBal_deltaLabel;

	/// Label for load balancing steps. 
	QLabel* loadBal_stepLabel;

	/// LineEdit for load balancing starttime.
	QLineEdit* loadBal_tLineEdit;

	/// LineEdit for load balancing timedelta.
	QLineEdit* loadBaL_deltaLineEdit;

	/// LineEdit for load balancing steps. 
	QLineEdit* loadBal_stepLineEdit;
          
	/// SolParallel_Data object storing the data. 
	SolParallel_Data* myDataPart;

	/// Configuration class for this SolParallel_GUI object. 
	TremoloGUIConfig* myConfig;

     
private slots:
     
};
 
#endif
