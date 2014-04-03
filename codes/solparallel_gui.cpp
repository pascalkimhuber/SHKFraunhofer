/***************************************************************************
                          solparallel_gui.cpp  -  description
                             -------------------
    begin                : Fri Jun 04 2004
			   Written by Daniel Bloemer, Tim Golla
    email                : bloemer@iam.uni-bonn.de, golla@ins.uni-bonn.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   Tremolo-X Project                                                     *
 *   All rights reserved                                                   *
 *   For details please see the LICENSE file                               *
 *                                                                         *
 ***************************************************************************/

// Standard-Includes 
#include <cfloat>


#include "solparallel_gui.h"
#include "solparallel_data.h"
#include "tremologuiconfig.h"
#include "qcomboboxwr.h"

// QT-Includes
#include <Q3ButtonGroup>
#include <QLabel>
#include <QLayout>
#include <QLCDNumber>
#include <QLineEdit>
#include <QRadioButton>
#include <QSpinBox>
#include <QValidator>
#include <Q3WidgetStack>
#include <Q3HBoxLayout>
#include <QGridLayout>
#include <Q3VBoxLayout>
#include <QToolTip>

SolParallel_GUI::SolParallel_GUI( SolParallel_Data* DataPart, TremoloGUIConfig* configIn, QWidget* parent, const char* name, Qt::WFlags f )
	:QWidget( parent, name, f )
{
	// Initialize data for this SolParallel_GUI object. 
	myDataPart = DataPart;
	// Update data if a hasChanged() signal is emitted.
	connect( myDataPart, SIGNAL( hasChanged() ), this, SLOT( updateGUI_data() ) );

	// Initialize configuration for this SolParallel_GUI object.
	myConfig = configIn;
	// Update data and labelsif a newScaling() signal is emitted. 
	connect( myConfig, SIGNAL( newScaling() ), this, SLOT( updateGUI_all() ) );
    
	// Initialize central layout. 
	centralLayout = new Q3VBoxLayout(this);
     
	// Set widgets and configuration for linked cell r_cut dialog (the text in the label is set in updateGUI_all). 
	Q3HBoxLayout* linkedCellR_CutLayout=new Q3HBoxLayout(centralLayout);
	linkedCellR_CutLabel = new QLabel( this, "linkedCellR_CutLabel" );
	linkedCellR_CutLayout->addWidget(linkedCellR_CutLabel);
	linkedCellR_CutLineEdit = new QLineEdit( QString::number( myDataPart->isLinkedCellR_Cut(), 'g' ), this, "linkedCellR_CutLineEdit" );
	linkedCellR_CutLayout->addWidget(linkedCellR_CutLineEdit);
	linkedCellR_CutLayout->addStretch();
	linkedCellR_CutLabel->setToolTip("Width of the linked cell grid.  "
									 "Changing this does not directly affect\n"
									 "the cutoff radius of short-range potentials, they will be computed\n"
									 "across several cells if necessary.");
	linkedCellR_CutLineEdit->setFixedWidth(200);
		     
	// Set validator for range checking for the linked cell r_cut LineEdit.
	QDoubleValidator* linkedCellR_CutLineEditValidator = new QDoubleValidator( DBL_MIN, DBL_MAX, 8, linkedCellR_CutLineEdit );
	linkedCellR_CutLineEdit->setValidator( linkedCellR_CutLineEditValidator );
	connect( linkedCellR_CutLineEdit, SIGNAL( returnPressed() ), this, SLOT( lCR_CLErp() ) );
	connect( linkedCellR_CutLineEdit, SIGNAL( lostFocus() ), this, SLOT( lCR_CLErp() ) );
     
	// Set widgets and configuration for coulomb force constant's dialog (the text for the label is set in updateGUI_all).
	Q3HBoxLayout* epsilon0invLayout=new Q3HBoxLayout(centralLayout);
	epsilon0invLabel = new QLabel( this, "epsilon0invLabel" );
	epsilon0invLayout->addWidget(epsilon0invLabel);
	epsilon0invLineEdit = new QLineEdit( QString::number( myDataPart->isEpsilon0inv(), 'g' ), this, "epsilon0invLineEdit" );
	epsilon0invLayout->addWidget(epsilon0invLineEdit);
	epsilon0invLayout->addStretch();
	epsilon0invLabel->setToolTip("Force constant before the Coulomb or gravitational force, in reduced units:\n"
								 "   epsilon0inv = 1/(4 pi epsilon0)  with epsilon0 = 8.854e-12 A^2 s^4 / kg m^3\n"
								 "For gravitation, set to -6.67428e-11 m^3 / kg s^2 and assign charges with mass converted with C/kg in SI units.");
	epsilon0invLineEdit->setFixedWidth(200);

	// Set validator for range checking for the coulomb force constant's LineEdit.
	QDoubleValidator* epsilon0invLineEditValidator = new QDoubleValidator( -DBL_MIN, DBL_MAX, 8, epsilon0invLineEdit );
	epsilon0invLineEdit->setValidator( epsilon0invLineEditValidator );
	connect( epsilon0invLineEdit, SIGNAL( returnPressed() ), this, SLOT( epsilon0inv_rp() ) );
	connect( epsilon0invLineEdit, SIGNAL( lostFocus() ), this, SLOT( epsilon0inv_rp() ) );

	// Add BoxLayout for the longrange solver's section. 
	Q3HBoxLayout* longRangeAndAdditionalValuesLayout=new Q3HBoxLayout(centralLayout);
     
	// Initialize buttongroup for the longrange solvers. 
	LongrangeButtonGroup = new Q3ButtonGroup( this, "LongrangeButtonGroup" );
	LongrangeButtonGroupLayout = new QGridLayout( LongrangeButtonGroup, 1, 1, 20, 5, "LongrangeButtonGroupLayout" );

	// Initialize radio buttons for all longrange solvers. 
	OFFRadioButton = new QRadioButton( LongrangeButtonGroup, "OFFRadioButton" );
	LongrangeButtonGroupLayout->addWidget( OFFRadioButton, 0, 0 );
	N_2RadioButton = new QRadioButton( LongrangeButtonGroup, "N_2RadioButton" );
	LongrangeButtonGroupLayout->addWidget( N_2RadioButton, 1, 0 );
	N_2SplineRadioButton = new QRadioButton( LongrangeButtonGroup, "N_2SplineRadioButton" );
	LongrangeButtonGroupLayout->addWidget( N_2SplineRadioButton, 2, 0 );
	SPMERadioButton = new QRadioButton( LongrangeButtonGroup, "SPMERadioButton" );
	LongrangeButtonGroupLayout->addWidget( SPMERadioButton, 3, 0 );
	FMMRadioButton = new QRadioButton( LongrangeButtonGroup, "FMMRadioButton" );
	LongrangeButtonGroupLayout->addWidget( FMMRadioButton, 4, 0 );
	FCS_DirectRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_DirectRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_DirectRadioButton, 5, 0);
	FCS_EwaldRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_EwaldRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_EwaldRadioButton, 6, 0);
	FCS_FMMRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_FMMRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_FMMRadioButton, 7, 0);
	FCS_PEPCRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_PEPCRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_PEPCRadioButton, 8, 0);
	FCS_PP3MGRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_PP3MGRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_PP3MGRadioButton, 9, 0);
	FCS_VMGRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_VMGRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_VMGRadioButton, 10, 0);
	FCS_P2NFFTRadioButton = new QRadioButton(LongrangeButtonGroup, "FCS_P2NFFTRadioButton");
	LongrangeButtonGroupLayout->addWidget(FCS_P2NFFTRadioButton, 11, 0);

	// Set the current selected longrange solver in the data class if a longrange solver's radio button has been clicked. 
	connect( LongrangeButtonGroup, SIGNAL( clicked ( int ) ), myDataPart, SLOT( setLongrangeAlgo( int ) ) );
	// Add the buttongroup for the longrange solvers to the BoxLayout. 
	longRangeAndAdditionalValuesLayout->addWidget(LongrangeButtonGroup);

	// Create buttongroup and layout for the parameters of the longrange solvers. 
	DataButtonGroup = new Q3ButtonGroup( this, "DataButtonGroup" );
	DataButtonGroup->setFixedWidth(400);
	DataButtonGroupLayout = new QGridLayout( DataButtonGroup, 1, 1, 20, 5, "DataButtonGroupLayout" );

	// Set spacing and alignment for DataButtonGroupLayout. 
	DataButtonGroupLayout->setColumnMinimumWidth(0, 200);     
	DataButtonGroupLayout->setAlignment(Qt::AlignTop);
     
	// Set Poisson solver section.
	poisson_solverLabel = new QLabel( DataButtonGroup, "poisson_solverLabel" );
	DataButtonGroupLayout->addWidget( poisson_solverLabel, 0, 0 );
	poisson_solverComboBox = new QComboBoxWR( DataButtonGroup, "poisson_solverComboBox" );
	poisson_solverComboBox->addItems( *myDataPart->isPoisson_solverList() );
	poisson_solverComboBox->setCurrentText( myDataPart->isPoisson_solver() );
	connect( poisson_solverComboBox, SIGNAL( activated ( const QString & ) ), myDataPart, SLOT( setPoisson_solver( const QString & ) ) );
	poisson_solverComboBox->setEnabled( false );
	DataButtonGroupLayout->addWidget( poisson_solverComboBox, 0, 1 );

	// Set r_cut parameter widgets.
	r_curLabel = new QLabel( DataButtonGroup, "r_curLabel" );
	DataButtonGroupLayout->addWidget( r_curLabel, 1, 0 );
	r_cutLineEdit = new QLineEdit( QString::number( myDataPart->isR_Cut(), 'g' ), DataButtonGroup, "r_cutLineEdit" );
	QDoubleValidator* r_cutLineEditEditValidator = new QDoubleValidator( DBL_MIN, DBL_MAX, 8, r_cutLineEdit );
	r_cutLineEdit->setValidator( r_cutLineEditEditValidator );
	connect( r_cutLineEdit, SIGNAL( returnPressed() ), this, SLOT( R_CLErp() ) );
	connect( r_cutLineEdit, SIGNAL( lostFocus() ), this, SLOT( R_CLErp() ) );
	r_cutLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( r_cutLineEdit, 1, 1 );

	// Set r_l parameter widgets. 
	r_lLabel = new QLabel( DataButtonGroup, "r_lLabel" );
	DataButtonGroupLayout->addWidget( r_lLabel, 2, 0 );
	r_lLineEdit = new QLineEdit( QString::number( myDataPart->isR_L(), 'g' ), DataButtonGroup, "r_lLineEdit" );
	QDoubleValidator* r_lLineEditValidator = new QDoubleValidator( DBL_MIN, DBL_MAX, 8, r_lLineEdit );
	r_lLineEdit->setValidator( r_lLineEditValidator );
	connect( r_lLineEdit, SIGNAL( returnPressed() ), this, SLOT( R_LLErp() ) );
	connect( r_lLineEdit, SIGNAL( lostFocus() ), this, SLOT( R_LLErp() ) );
	r_lLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( r_lLineEdit, 2, 1 );

	// Set splitting coefficient parameter widgets.
	splittingCoeffizientGLabel = new QLabel( DataButtonGroup, "splittingCoeffizientGLabel" );
	DataButtonGroupLayout->addWidget( splittingCoeffizientGLabel, 3, 0 );
	splittingCoeffizientGLineEdit = new QLineEdit( QString::number( myDataPart->isSplittingCoefficientG(), 'g' ), DataButtonGroup, "splittingCoeffizientGLineEdit" );
	QDoubleValidator* splittingCoeffizientGLineEditValidator = new QDoubleValidator( 0, 1, 8, splittingCoeffizientGLineEdit );
	splittingCoeffizientGLineEdit->setValidator( splittingCoeffizientGLineEditValidator );
	connect( splittingCoeffizientGLineEdit, SIGNAL( returnPressed() ), this, SLOT( sCGLErp() ) );
	connect( splittingCoeffizientGLineEdit, SIGNAL( lostFocus() ), this, SLOT( sCGLErp() ) );
	splittingCoeffizientGLineEdit->setEnabled( false );

	// Set MAP parameter widgets.
	DataButtonGroupLayout->addWidget( splittingCoeffizientGLineEdit, 3, 1 );
	MapLabel = new QLabel( DataButtonGroup, "MapLabel" );
	DataButtonGroupLayout->addWidget( MapLabel, 4, 0 );
	MapLineEdit = new QLineEdit( QString::number( myDataPart->isMAP(), 'g' ), DataButtonGroup, "MapLineEdit" );
	QDoubleValidator* MapLineEditValidator = new QDoubleValidator( 0, 1, 8, MapLineEdit );
	MapLineEdit->setValidator( MapLineEditValidator );
	connect( MapLineEdit, SIGNAL( returnPressed() ), this, SLOT( mapLErp() ) );
	connect( MapLineEdit, SIGNAL( lostFocus() ), this, SLOT( mapLErp() ) );
	MapLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( MapLineEdit, 4, 1 );

	// Set cellratio parameter widgets. 
	cellratioLabel = new QLabel( DataButtonGroup, "cellratioLabel" );
	DataButtonGroupLayout->addWidget( cellratioLabel, 5, 0 );
	cellratioLineEdit = new QLineEdit( QString::number( myDataPart->isCellratio() ), DataButtonGroup, "cellratioLineEdit" );
	QIntValidator* cellratioLineEditValidator = new QIntValidator( 0, INT_MAX, cellratioLineEdit );
	cellratioLineEdit->setValidator( cellratioLineEditValidator );
	connect( cellratioLineEdit, SIGNAL( returnPressed() ), this, SLOT( cellratioLErp() ) );
	connect( cellratioLineEdit, SIGNAL( lostFocus() ), this, SLOT( cellratioLErp() ) );
	cellratioLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( cellratioLineEdit, 5, 1 );

	// Set interpolation degree parameter widgets. 
	interpolationLabel = new QLabel( DataButtonGroup, "interpolationLabel" );
	DataButtonGroupLayout->addWidget( interpolationLabel, 6, 0 );
	interpolationLineEdit = new QLineEdit( QString::number( myDataPart->isInterpolationDegree() ), DataButtonGroup, "interpolationLineEdit" );
	QIntValidator* interpolationLineEditValidator = new QIntValidator( 0, INT_MAX, interpolationLineEdit );
	interpolationLineEdit->setValidator( interpolationLineEditValidator );
	connect( interpolationLineEdit, SIGNAL( returnPressed() ), this, SLOT( interpolationLErp() ) );
	connect( interpolationLineEdit, SIGNAL( lostFocus() ), this, SLOT( interpolationLErp() ) );
	interpolationLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( interpolationLineEdit, 6, 1 );

	// Set maximum tree level parameter widgets. 
	maxTreeLabel = new QLabel( DataButtonGroup, "maxTreeLabel" );
	DataButtonGroupLayout->addWidget( maxTreeLabel, 7, 0 );
	maxTreeLineEdit = new QLineEdit( QString::number( myDataPart->isMaxTreeLevel() ), DataButtonGroup, "maxTreeLineEdit" );
	QIntValidator* maxTreeLineEditValidator = new QIntValidator( 0, INT_MAX, maxTreeLineEdit );
	maxTreeLineEdit->setValidator( maxTreeLineEditValidator );
	connect( maxTreeLineEdit, SIGNAL( returnPressed() ), this, SLOT( maxTreeLErp() ) );
	connect( maxTreeLineEdit, SIGNAL( lostFocus() ), this, SLOT( maxTreeLErp() ) );
	maxTreeLineEdit->setEnabled( false );
	DataButtonGroupLayout->addWidget( maxTreeLineEdit, 7, 1 );

	// Set splitting coefficient alpha widgets. 
	fcsSplittingCoefficientAlphaLabel = new QLabel(DataButtonGroup, "fcsSplittingCoefficientAlphaLabel");
	DataButtonGroupLayout->addWidget(fcsSplittingCoefficientAlphaLabel, 8, 0);
	fcsSplittingCoefficientAlphaLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_splittingCoefficientAlpha()), DataButtonGroup, "fcsSplittingCoefficientAlphaLieEdit");
	QDoubleValidator* fcsSplittingCoefficientAlphaLineEditValidator = new QDoubleValidator(DBL_MIN, DBL_MAX, 8, fcsSplittingCoefficientAlphaLineEdit);
	fcsSplittingCoefficientAlphaLineEdit->setValidator(fcsSplittingCoefficientAlphaLineEditValidator);
	connect(fcsSplittingCoefficientAlphaLineEdit, SIGNAL(returnPressed()), this, SLOT(fcssCAlphaLErp()));
	connect(fcsSplittingCoefficientAlphaLineEdit, SIGNAL(lostFocus()), this, SLOT(fcssCAlphaLErp()));
	fcsSplittingCoefficientAlphaLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsSplittingCoefficientAlphaLineEdit, 8, 1);

	// Set tolerance widges. 
	fcsToleranceLabel = new QLabel(DataButtonGroup, "fcsToleranceLabel");
	DataButtonGroupLayout->addWidget(fcsToleranceLabel, 9, 0);
	fcsToleranceLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_tolerance(), 'g'), DataButtonGroup, "fcsToleranceLineEdit");
	

	// /// Label for tolerance (ScaFaCoS). 
	// QLabel* fcsToleranceLabel;

	// /// LineEdit for tolerance (ScaFaCoS). 
	// QLineEdit* fcsToleranceLineEdit;

	// /// Label for tolerance_type (ScaFaCoS). 
	// QLabel* fcsToleranceTypeLabel;
	
	// /// LineEdit for tolerance_type (ScaFaCoS). 
	// QLineEdit* fcsToleranceTypeLineEdit;

	// /// Label for periodic_images (ScaFaCoS).
	// QLabel* fcsPeriodicImagesLabel;
	
	// /// LineEdit for periodic_images_x (ScaFaCoS). 
	// QLineEdit* fcsPeriodicImagesXLineEdit;

	// /// LineEdit for periodic_images_y (ScaFaCoS). 
	// QLineEdit* fcsPeriodicImagesYLineEdit;

	// /// LineEdit for periodic_images_z (ScaFaCoS). 
	// QLineEdit* fcsPeriodicImagesZLineEdit;

	// /// Label for kmax (ScaFaCoS). 
	// QLabel* fcsKmaxLabel;
	
	// /// LineEdit for kmax (ScaFaCoS). 
	// QLineEdit* fcsKmaxLineEdit;

	// /// Label for maxkmax (ScaFaCoS). 
	// QLabel* fcsMaxkmaxLabel;

	// /// LineEdit for maxkmax (ScaFaCoS).
	// QLineEdit* fcsMaxkmaxLineEdit;

	// /// Label for balanceload (ScaFaCoS). 
	// QLabel* fcsBalanceloadLabel;

	// /// LineEdit for balanceload (ScaFaCoS). 
	// QLineEdit* fcsBalanceloadLineEdit;

	// /// Label for dipole_correction (ScaFaCoS). 
	// QLabel* fcsDipoleCorrectionLabel;
	
	// /// LineEdit for dipole_correction (ScaFaCoS). 
	// QLineEdit* fcsDipoleCorrectionLineEdit;

	// /// Label for maxdepth (ScaFaCoS). 
	// QLabel* fcsMaxdepthLabel;

	// /// LineEdit for maxdepth (ScaFaCoS). 
	// QLineEdit* fcsMaxdepthLineEdit;

	// /// Label for potential (ScaFaCoS). 
	// QLabel* fcsPotentialLabel;

	// /// LineEdit for potential (ScaFaCoS). 
	// QLineEdit* fcsPotentialLineEdit;

	// /// Label for radius (ScaFaCoS). 
	// QLabel* fcsRadiusLabel;
	
	// /// LineEdit for radius (ScaFaCoS). 
	// QLineEdit* fcsRadiusLineEdit;

	// /// Label for unroll_limit (ScaFaCoS). 
	// QLabel* fcsUnrollLimitLabel;

	// /// LineEdit for unroll_limit (ScaFaCoS). 
	// QLineEdit* fcsUnrollLimitLineEdit;

	// /// Label for degree (ScaFaCoS). 
	// QLabel* fcsDegreeLabel;

	// /// LineEdit for degree (ScaFaCoS). 
	// QLineEdit* fcsDegreeLineEdit;

	// /// Label for ghosts (ScaFaCoS). 
	// QLabel* fcsGhostsLabel;

	// /// LineEdit for ghosts (ScaFaCoS). 
	// QLineEdit* fcsGhostsLineEdit;

	// /// Label for gridsize (ScaFaCoS). 
	// QLabel* fcsGridsizeLabel;

	// /// LineEdit for gridsize_x (ScaFaCoS). 
	// QLineEdit* fcsGridsizeXLineEdit;

	// /// LineEdit for gridsize_y (ScaFaCoS). 
	// QLineEdit* fcsGridsizeYLineEdit;

	// /// LineEdit for gridsize_z (ScaFaCoS). 
	// QLineEdit* fcsGridsizeZLineEdit;

	// /// Label for max_iterations (ScaFaCoS). 
	// QLabel* fcsMaxIterationsLabel;

	// /// LineEdit for max_iterations (ScaFaCoS). 
	// QLineEdit* fcsMaxInterationsLineEdit;

	// /// Label for debuglevel (ScaFaCoS). 
	// QLabel* fcsDebuglevelLabel;

	// /// LineEdit for debuglevel (ScaFaCoS). 
	// QLineEdit* fcsDebuglevelLineEdit;

	// /// Label for epsilon (ScaFaCoS). 
	// QLabel* fcsEpsilonLabel;

	// /// LineEdit for epsilon (ScaFaCoS). 
	// QLineEdit* fcsEpsilonLineEdit;

	// /// Label for load_balancing (ScaFaCoS). 
	// QLabel* fcsLoadBalancingLabel;

	// /// LineEdit for load_balancing (ScaFaCoS). 
	// QLineEdit* fcsLoadBalancingLineEdit;

	// /// Label for npm (ScaFaCoS). 
	// QLabel* fcsNpmLabel;

	// /// LineEdit for npm (ScaFaCoS). 
	// QLineEdit* fcsNpmLineEdit;

	// /// Label for num_walk_threads (ScaFaCoS). 
	// QLabel* fcsNumWalkThreadsLabel;

	// /// LineEdit for num_walk_threads (ScaFaCoS). 
	// QLineEdit* fcsNumWalkThreadsLineEdit;

	// /// Label for theta (ScaFaCoS). 
	// QLabel* fcsThetaLabel;

	// /// LineEdit for theta (ScaFaCoS). 
	// QLineEdit* fcsThetaLineEdit;

	// /// Label for cycle_type (ScaFaCoS). 
	// QLabel* fcsCycleTypeLabel;

	// /// LineEdit for cycle_type (ScaFaCoS). 
	// QLineEdit* fcsCycleTypeLineEdit;

	// /// Label for discretization_order (ScaFaCoS). 
	// QLabel* fcsDiscretizationOrderLabel;

	// /// LineEdit for discretization_order (ScaFaCoS). 
	// QLineEdit* fcsDiscretizationOrderLineEdit;

	// /// Label for interpolation_order (ScaFaCoS). 
	// QLabel* fcsInterpolationOrderLabel;

	// /// LineEdit for interpolation_order (ScaFaCoS). 
	// QLineEdit* fcsInterpolationOrderLineEdit;

	// /// Label for max_level (ScaFaCoS). 
	// QLabel* fcsMaxLevelLabel;

	// /// LineEdit for max_level (ScaFaCoS). 
	// QLineEdit* fcsMaxLevelLineEdit;

	// /// Label for near_field_cells (ScaFaCoS). 
	// QLabel* fcsNearFieldCellsLabel;

	// /// LineEdit for near_field_cells (ScaFaCoS). 
	// QLineEdit* fcsNearFieldCellsLineEdit;

	// /// Label for precision (ScaFaCoS). 
	// QLabel* fcsPrecisionLabel;
	
	// /// LineEdit for precision (ScaFaCoS). 
	// QLineEdit* fcsPrecisionLineEdit;

	// /// Label for smoothing_steps (ScaFaCoS). 
	// QLabel* fcsSmoothingStepsLabel;

	// /// LineEdit for smoothing_steps (ScaFaCoS). 
	// QLineEdit* fcsSmoothingStepsLineEdit;

	// /// Label for epsI (ScaFaCoS). 
	// QLabel* fcsEpsILabel;

	// /// LineEdit for epsI (ScaFaCoS). 
	// QLineEdit* fcsEpsILineEdit;

	// /// Label for m (ScaFaCoS). 
	// QLabel* fcsMLabel;

	// /// LineEdit for m (ScaFaCoS). 
	// QLineEdit* fcsMLineEdit;

	// /// Label for p (ScaFaCoS). 
	// QLabel* fcsPLabel;

	// /// LineEdit for p (ScaFaCoS). 
	// QLineEdit* fcsPLineEdit;

	// /// Label for oversampled_gridsize (ScaFaCoS). 
	// QLabel* fcsOversampledGridsizeLabel;

	// /// LineEdit for oversampled_gridsize_x (ScaFaCoS). 
	// QLineEdit* fcsOversampledGridsizeXLineEdit;

	// /// LineEdit for oversampled_gridsize_y (ScaFaCoS). 
	// QLineEdit* fcsOversampledGridsizeYLineEdit;

	// /// LineEdit for oversampled_gridsize_z (ScaFaCoS). 
	// QLineEdit* fcsOversampledGridsizeZLineEdit;

     
	// Configure visibility and enabling of poisson solver widgets.
	connect( myDataPart, SIGNAL( enablepoisson_solver( bool ) ), poisson_solverComboBox, SLOT( setEnabled( bool ) ) );
	connect(myDataPart, SIGNAL(enablepoisson_solver(bool)), poisson_solverComboBox, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enablepoisson_solver(bool)), poisson_solverLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of r_cut widgets. 
	connect( myDataPart, SIGNAL( enableR_Cut( bool ) ), r_cutLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableR_Cut( bool ) ), r_cutLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableR_Cut( bool ) ), r_curLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of r_l widgets. 
	connect( myDataPart, SIGNAL( enableR_L( bool ) ), r_lLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableR_L( bool ) ), r_lLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableR_L( bool ) ), r_lLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of splitting coefficient G widgets. 
	connect( myDataPart, SIGNAL( enableSplittingCoefficientG( bool ) ), splittingCoeffizientGLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableSplittingCoefficientG( bool ) ), splittingCoeffizientGLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableSplittingCoefficientG( bool ) ), splittingCoeffizientGLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of MAP widgets. 
	connect( myDataPart, SIGNAL( enableMAP( bool ) ), MapLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableMAP( bool ) ), MapLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableMAP( bool ) ), MapLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of maximum tree level widgets. 
	connect( myDataPart, SIGNAL( enableMaxTreeLevel( bool ) ), maxTreeLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableMaxTreeLevel( bool ) ), maxTreeLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableMaxTreeLevel( bool ) ), maxTreeLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of interpolation degree widgets. 
	connect( myDataPart, SIGNAL( enableInterpolationDegree( bool ) ), interpolationLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableInterpolationDegree( bool ) ), interpolationLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableInterpolationDegree( bool ) ), interpolationLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of cellration widgets. 
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLabel, SLOT( setVisible( bool ) ) );
     
	// Add buttongroup for the parameters to the longrange solver's layout. 
	longRangeAndAdditionalValuesLayout->addWidget(DataButtonGroup);
	longRangeAndAdditionalValuesLayout->addStretch();

	// Initialize layout and buttongroup for the parallelization section. 
	Q3HBoxLayout* parallelizerButtonGroupHBoxLayout=new Q3HBoxLayout(centralLayout);
	parallelizerButtonGroup = new Q3ButtonGroup( this, "parallelizerButtonGroup" ); 

	// Add parallelization buttongroup to the parallelization layout. 
	parallelizerButtonGroupHBoxLayout->addWidget(parallelizerButtonGroup);
	parallelizerButtonGroupHBoxLayout->addStretch();

	// Add GritLayout to parallelization ButtonGroup. 
	parallelizerButtonGroupLayout = new QGridLayout( parallelizerButtonGroup, 1, 1, 20, 5, "parallelizerButtonGroupLayout" );

	// Set labels for load balancing. 
	loadBalancingfunctionLabel = new QLabel( parallelizerButtonGroup, "loadBalancingfunctionLabel" );
	parallelizerButtonGroupLayout->addWidget( loadBalancingfunctionLabel, 0, 0 );

	cutsXLabel = new QLabel( parallelizerButtonGroup, "cutsXLabel" );
	parallelizerButtonGroupLayout->addWidget( cutsXLabel, 0, 1 );

	cutsYLabel = new QLabel( parallelizerButtonGroup, "cutsYLabel" );
	parallelizerButtonGroupLayout->addWidget( cutsYLabel, 0, 2 );

	cutsZLabel = new QLabel( parallelizerButtonGroup, "cutsZLabel" );
	parallelizerButtonGroupLayout->addWidget( cutsZLabel, 0, 3 );

	usedProcsLabel = new QLabel( parallelizerButtonGroup, "usedProcsLabel" );
	parallelizerButtonGroupLayout->addWidget( usedProcsLabel, 0, 4 );

	// Set ComboBoxes and SpinBoxes for load balancing. 
	loadBalancingfunctionComboBox = new QComboBoxWR( parallelizerButtonGroup, "loadBalancingfunctionComboBox" );
	loadBalancingfunctionComboBox->addItems( *myDataPart->isLoadBalancingFunctionList() );
	loadBalancingfunctionComboBox->setCurrentText( myDataPart->isLoadBalancingFunction() );
	connect( loadBalancingfunctionComboBox, SIGNAL( activated ( const QString & ) ), myDataPart, SLOT( setLoadBalancingFunction( const QString & ) ) );
	parallelizerButtonGroupLayout->addWidget( loadBalancingfunctionComboBox, 1, 0 );

	cutsXSpinBox = new QSpinBox( parallelizerButtonGroup, "cutsXSpinBox" );
	cutsXSpinBox->setMinValue( 1 );
	cutsXSpinBox->setMaxValue( INT_MAX );
	cutsXSpinBox->setButtonSymbols( QSpinBox::PlusMinus );
	cutsXSpinBox->setValue( myDataPart->isCutsX() );
	connect( cutsXSpinBox, SIGNAL( valueChanged( int ) ), myDataPart, SLOT( setCutsX( int  ) ) );
	parallelizerButtonGroupLayout->addWidget( cutsXSpinBox, 1, 1 );

	cutsYSpinBox = new QSpinBox( parallelizerButtonGroup, "cutsYSpinBox" );
	cutsYSpinBox->setMinValue( 1 );
	cutsYSpinBox->setMaxValue( INT_MAX );
	cutsYSpinBox->setButtonSymbols( QSpinBox::PlusMinus );
	cutsYSpinBox->setValue( myDataPart->isCutsY() );
	connect( cutsYSpinBox, SIGNAL( valueChanged( int ) ), myDataPart, SLOT( setCutsY( int  ) ) );
	parallelizerButtonGroupLayout->addWidget( cutsYSpinBox, 1, 2 );

	cutsZSpinBox = new QSpinBox( parallelizerButtonGroup, "cutsZSpinBox" );
	cutsZSpinBox->setMinValue( 1 );
	cutsZSpinBox->setMaxValue( INT_MAX );
	cutsZSpinBox->setButtonSymbols( QSpinBox::PlusMinus );
	cutsZSpinBox->setValue( myDataPart->isCutsZ() );
	connect( cutsZSpinBox, SIGNAL( valueChanged( int ) ), myDataPart, SLOT( setCutsZ( int  ) ) );
	parallelizerButtonGroupLayout->addWidget( cutsZSpinBox, 1, 3 );

	// Add SpinBox and value label to WidgetStack. By this either the total number of processes is displayed or the number of processes can be set. 
	usedProcsWidgetStack = new Q3WidgetStack( parallelizerButtonGroup, "usedProcsWidgetStack" );
	usedProcsValueLabel = new QLabel( usedProcsWidgetStack, "usedProcsValueLabel" );
	usedProcsValueLabel->setAlignment(Qt::AlignRight);
	usedProcsWidgetStack->addWidget( usedProcsValueLabel, 0 );

	// Configure SpinBox for the number of processes. 
	usedProcsSpinBox = new QSpinBox( usedProcsWidgetStack, "usedProcsSpinBox" );
	usedProcsSpinBox->setMinValue( 1 );
	usedProcsSpinBox->setMaxValue( INT_MAX );
	usedProcsSpinBox->setButtonSymbols( QSpinBox::PlusMinus );
	usedProcsSpinBox->setValue( myDataPart->isProcsAtAll() );
	connect( usedProcsSpinBox, SIGNAL( valueChanged( int ) ), myDataPart, SLOT( setProcsAtAll( int  ) ) );
	usedProcsWidgetStack->addWidget( usedProcsSpinBox, 1 );
	usedProcsWidgetStack->raiseWidget( 0 );

	// Add Labels for starttime, timedelta and steps for load balancing. 
	parallelizerButtonGroupLayout->addWidget( usedProcsWidgetStack, 1, 4 );
	loadBal_tLabel = new QLabel( parallelizerButtonGroup, "loadBal_tLabel" );
	parallelizerButtonGroupLayout->addWidget( loadBal_tLabel, 2, 0 );
	loadBal_deltaLabel = new QLabel( parallelizerButtonGroup, "loadBal_deltaLabel" );
	parallelizerButtonGroupLayout->addWidget( loadBal_deltaLabel, 3, 0 );
	loadBal_stepLabel = new QLabel( parallelizerButtonGroup, "loadBal_stepLabel" );
	parallelizerButtonGroupLayout->addWidget( loadBal_stepLabel, 4, 0 );

	// Add LineEdits for starttime, timedelta and steps for load balancing. 
	loadBal_tLineEdit = new QLineEdit( QString::number( myDataPart->isLoadbal_t(), 'g' ), parallelizerButtonGroup, "loadBal_tLineEdit" );
	QDoubleValidator* loadBal_tLineEditValidator = new QDoubleValidator( 0, DBL_MAX, 8, loadBal_tLineEdit );
	loadBal_tLineEdit->setValidator( loadBal_tLineEditValidator );
	connect( loadBal_tLineEdit, SIGNAL( returnPressed() ), this, SLOT( lb_tLErp() ) );
	connect( loadBal_tLineEdit, SIGNAL( lostFocus() ), this, SLOT( lb_tLErp() ) );
	loadBal_tLineEdit->setEnabled( false );
	parallelizerButtonGroupLayout->addMultiCellWidget( loadBal_tLineEdit, 2, 2, 1, 4 );

	loadBaL_deltaLineEdit = new QLineEdit( QString::number( myDataPart->isLoadbal_delta(), 'g' ), parallelizerButtonGroup, "loadBaL_deltaLineEdit" );
	QDoubleValidator* loadBaL_deltaLineEditValidator = new QDoubleValidator( DBL_MIN, DBL_MAX, 8, loadBaL_deltaLineEdit );
	loadBaL_deltaLineEdit->setValidator( loadBaL_deltaLineEditValidator );
	connect( loadBaL_deltaLineEdit, SIGNAL( returnPressed() ), this, SLOT( lb_deltaLErp() ) );
	connect( loadBaL_deltaLineEdit, SIGNAL( lostFocus() ), this, SLOT( lb_deltaLErp() ) );
	loadBaL_deltaLineEdit->setEnabled( false );
	parallelizerButtonGroupLayout->addMultiCellWidget( loadBaL_deltaLineEdit, 3, 3, 1, 4 );

	loadBal_stepLineEdit = new QLineEdit( QString::number( myDataPart->isLoadbal_step() ), parallelizerButtonGroup, "loadBal_stepLineEdit" );
	QIntValidator* loadBal_stepLineEditValidator = new QIntValidator( 1, INT_MAX, loadBal_stepLineEdit );
	loadBal_stepLineEdit->setValidator( loadBal_stepLineEditValidator );
	connect( loadBal_stepLineEdit, SIGNAL( returnPressed() ), this, SLOT( lb_stepLErp() ) );
	connect( loadBal_stepLineEdit, SIGNAL( lostFocus() ), this, SLOT( lb_stepLErp() ) );
	loadBal_stepLineEdit->setEnabled( false );
	parallelizerButtonGroupLayout->addMultiCellWidget( loadBal_stepLineEdit, 4, 4, 1, 4 );

	// Connect enable-signals to enable-slots for the LineEdits of the load balancing LineEdits.
	connect( myDataPart, SIGNAL( enableLoadBalValues_dyn( bool ) ), loadBal_tLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableLoadBalValues_dyn( bool ) ), loadBaL_deltaLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableLoadBalValues_opt( bool ) ), loadBal_stepLineEdit, SLOT( setEnabled( bool ) ) );
     
	centralLayout->addStretch();

	// Call update-methods to update displayed data. 
	updateGUI_all();
	languageChange(); 
}

SolParallel_GUI::~SolParallel_GUI()
{
}


void SolParallel_GUI::languageChange()
{
	// Set titles in right language. 
    setCaption( tr( "SolParallel_GUI" ) );
    LongrangeButtonGroup->setTitle( tr( "Longrange Algorithms:" ) );

	// Set text for the longrange solver's labels. 
    OFFRadioButton->setText( tr( "Off" ) );
    QString N2String("N");
    N2String.append(TremoloGUIConfig::numberToUnicodeExponent( 2));
//     N_2RadioButton->setText( tr( "N^2" ) );
    N_2RadioButton->setText( N2String );
    N2String.append(" Spline");
//     N_2SplineRadioButton->setText( tr( "N^2 Spline" ) );
    N_2SplineRadioButton->setText( N2String );
    SPMERadioButton->setText( tr( "SPME" ) );
    FMMRadioButton->setText( tr( "FMM" ) );
	FCS_DirectRadioButton->setText(tr("Direct summation (ScaFaCoS)"));
	FCS_EwaldRadioButton->setText(tr("Ewald (ScaFaCoS)"));
	FCS_FMMRadioButton->setText(tr("FMM (ScaFaCoS)"));
	FCS_PEPCRadioButton->setText(tr("PEPC (ScaFaCoS)"));
	FCS_PP3MGRadioButton->setText(tr("PP3MG (ScaFaCoS)"));
	FCS_VMGRadioButton->setText(tr("VMG (ScaFaCoS)"));
	FCS_P2NFFTRadioButton->setText(tr("P2NFFT (ScaFaCoS)"));

	// Set text for the parameters labels. 
    DataButtonGroup->setTitle( tr( "Additional Values:" ) );
    poisson_solverLabel->setText( tr( "Poisson Solver:" ) );
    interpolationLabel->setText( tr( "Interpolation Degree:" ) );
    maxTreeLabel->setText( tr( "Maximum tree level:" ) );
    cellratioLabel->setText( tr( "Cellratio:" ) );
	fcsToleranceLabel->setText(tr("Tolerance:"));
	fcsToleranceTypeLabel->setText(tr("Tolerance type:"));
	fcsPeriodicImagesLabel->setText(tr("Number of periodic images per axis:"));
	fcsKmaxLabel->setText(tr("Number of wave vectors (kmax):"));
	fcsMaxkmaxLabel->setText(tr("Maximum number of wave vectors (maxkmax):"));
	fcsBalanceloadLabel->setText(tr("Load balancing routine:"));
	fcsDipoleCorrectionLabel->setText(tr("Dipole correction:"));
	fcsMaxdepthLabel->setText(tr("Maximum tree depth:"));
	fcsPotentialLabel->setText(tr("Potential:"));
	fcsRadiusLabel->setText(tr("Radius for cusp potential:"));
	fcsUnrollLimitLabel->setText(tr("Limit for unrolled functions:"));
	fcsDegreeLabel->setText(tr("Interpolation degree:"));
	fcsGhostsLabel->setText(tr("Number of ghosts cells:"));
	fcsGridsizeLabel->setText(tr("Number of grid cells per axis:"));
	fcsMaxIterationsLabel->setText(tr("Maximum number of iterations:"));
	fcsDebuglevelLabel->setText(tr("Debug level:"));
	fcsEpsilonLabel->setText(tr("epsilon:"));
	fcsNpmLabel->setText(tr("NPM:"));
	fcsNumWalkThreadsLabel->setText(tr("Number of walk threads per MPI rank:"));
	fcsThetaLabel->setText(tr("theta:"));
	fcsCycleTypeLabel->setText(tr("Multi-grid cycle type:"));
	fcsDiscretizationOrderLabel->setText(tr("Order of discretization scheme for PDE solver:"));
	fcsInterpolationOrderLabel->setText(tr("Degree of interpolation polynomial:"));
	fcsMaxLevelLabel->setText(tr("Maximum number of multi-grid levels:"));
	fcsNearFieldCellsLabel->setText(tr("Range of smearing of charges over near field cells:"));
	fcsPrecisionLabel->setText(tr("Threshold for residual for solver iteration:"));
	fcsSmoothingStepsLabel->setText(tr("Number of smoothing steps:"));
	fcsEpsILabel->setText(tr("epsI:"));
	fcsMLabel->setText(tr("m:"));
	fcsPLabel->setText(tr("p:"));
	fcsOversampledGridsizeLabel->setText(tr("Size of oversampled FFT grid:"));
	
	// Set text for the parallelization labels. 
    parallelizerButtonGroup->setTitle( tr( "Parallelization" ) );
	loadBalancingfunctionLabel->setText( tr( "LoadBalancing:" ) );
	cutsXLabel->setText( tr( "Cuts in X:" ) );
	cutsYLabel->setText( tr( "Cuts in Y:" ) );
	cutsZLabel->setText( tr( "Cuts in Z:" ) );
	usedProcsLabel->setText( tr( "Used Processes:" ) );
	loadBal_stepLabel->setText( tr( "LB step:" ) );
}

void SolParallel_GUI::updateGUI_all()
{
	// Update units in all labels. 
    linkedCellR_CutLabel->setText( tr( "Linked Cell r_cut (" + myDataPart->isLinkedCellR_CutUnit() + "):" ) );
    epsilon0invLabel->setText( tr( "Coulomb force constant  (" + myDataPart->isEpsilon0invUnit() + "):" ) );

    r_curLabel->setText( tr( "r_cut (" + myDataPart->isR_CutUnit() + "):" ) );
    r_lLabel->setText( tr( "r_l (" + myDataPart->isR_LUnit() + "):" ) );
    splittingCoeffizientGLabel->setText( tr( "Splitting Coefficient G (" + myDataPart->isSplittingCoefficientGUnit() + "):" ) );
    MapLabel->setText( tr( "MAP (" + myDataPart->isMAPUnit() + "):" ) );

	QString alpha;
	alpha = QChar(0xb1, 0x03);
	fcsSplittingCoefficientAlphaLabel->setText(tr("Splitting Coefficient " + alpha + " (" + myDataPart->isFCS_splittingCoefficientAlphaUnit() + "):"));

    loadBal_tLabel->setText( tr( "LB Starttime (" + myDataPart->isLoadbal_tUnit() + "):" ) );
    loadBal_deltaLabel->setText( tr( "LB Timedelta (" + myDataPart->isLoadbal_deltaUnit() + "):" ) );

	// Update data. 
    updateGUI_data();
}

void SolParallel_GUI::updateGUI_data()
{
	// Update LineEdit for "Linked Cell r_cut".
	linkedCellR_CutLineEdit->setText( QString::number( myDataPart->isLinkedCellR_Cut(), 'g' ) );

	// Update LineEdit for "Coulomb force constant". 
	epsilon0invLineEdit->setText( QString::number( myDataPart->isEpsilon0inv(), 'g' ) );

	// Update RadioButton indicating the current longrange solver. 
	LongrangeButtonGroup->setButton( myDataPart->isLongrangeAlgo() );

	// Update ComboBox for "Poisson solver". 
	poisson_solverComboBox->setCurrentText( myDataPart->isPoisson_solver() );

	// Update LineEdit for "r_cut" of the current solver.
	r_cutLineEdit->setText( QString::number( myDataPart->isR_Cut(), 'g' ) );

	// Update LineEdit for "r_l" of the current solver. 
	r_lLineEdit->setText( QString::number( myDataPart->isR_L(), 'g' ) );

	// Update LineEdit for "Splitting Coefficient G" of the current solver.
	splittingCoeffizientGLineEdit->setText( QString::number( myDataPart->isSplittingCoefficientG(), 'g' ) );

	// Update LineEdit for "MAP ()" of the current solver. 
	MapLineEdit->setText( QString::number( myDataPart->isMAP(), 'g' ) );

	// Update LineEdit for "Cellratio" of the current solver. 
	cellratioLineEdit->setText( QString::number( myDataPart->isCellratio() ) ); // Noch anpassen....

	// Update LineEdit for "Interpolation Degree" of the current solver.
	interpolationLineEdit->setText( QString::number( myDataPart->isInterpolationDegree() ) );

	// Update LineEdit for "Maximum tree level" of the current solver. 
	maxTreeLineEdit->setText( QString::number( myDataPart->isMaxTreeLevel() ) );

	// Update LineEdit for "Splitting coefficient alpha". 
	fcsSplittingCoefficientAlphaLineEdit->setText(QString::number(myDataPart->isFCS_splittingCoefficientAlpha(), 'g'));

	// Update LineEdit for "Tolerance". 
	fcsToleranceLineEdit->setText(QString::number(myDataPart->isFCS_tolerance(), 'g'));
	
	/// Update LineEdit for "Tolerance type". 
	fcsToleranceTypeLineEdit->setText(QString::number(myDataPart->isFCS_tolerance_type()));

	/// Update LineEdits for for "Number of periodic images per axis". 
	fcsPeriodicImagesXLineEdit->setText(QString::number(myDataPart->isFCS_periodic_images_x()));
	fcsPeriodicImagesYLineEdit->setText(QString::number(myDataPart->isFCS_periodic_images_y()));
	fcsPeriodicImagesZLineEdit->setText(QString::number(myDataPart->isFCS_periodic_images_z()));

	/// Update LineEdit for "Number of wave vectors (kmax)". 
	fcsKmaxLineEdit->setText(QString::number(myDataPart->isFCS_kmax()));
	
	/// Update LineEdit for "Maximum knumber of wave vectors (maxkmax)". 
	fcsMaxkmaxLineEdit->setText(QString::number(myDataPart->isFCS_maxkmax()));

	/// Update LineEdit for "Load balancing routine". 
	fcsBalanceloadLineEdit->setText(QString::number(myDataPart->isFCS_balanceload()));

	/// Update LineEdit for "Dipole correction". 
	fcsDipoleCorrectionLineEdit->setText(QString::number(myDataPart->isFCS_dipole_correction()));

	/// Update LineEdit for "Maximum tree depth". 
	fcsMaxdepthLineEdit->setText(QString::number(myDataPart->isFCS_maxdepth()));

	/// Update LineEdit for "Potential". 
	fcsPotentialLineEdit->setText(QString::number(myDataPart->isFCS_potential()));

	/// Update LineEdit for "Radius for cusp potential". 
	fcsRadiusLineEdit->setText(QString::number(myDataPart->isFCS_radius()));

	/// Update LineEdit for "Limit for unrolled functions". 
	fcsUnrollLimitLineEdit->setText(QString::number(myDataPart->isFCS_unroll_limit()));

	/// Update LineEdit for "Interpolation degree". 
	fcsDegreeLineEdit->setText(QString::number(myDataPart->isFCS_degree()));

	/// Update LineEdit for "Number of ghosts cells". 
	fcsGhostsLineEdit->setText(QString::number(myDataPart->isFCS_ghosts()));
 
	/// Update LineEdits for "Number of grid cells per axis". 
	fcsGridsizeXLineEdit->setText(QString::number(myDataPart->isFCS_gridsize_x()));
	fcsGridsizeYLineEdit->setText(QString::number(myDataPart->isFCS_gridsize_y()));
	fcsGridsizeZLineEdit->setText(QString::number(myDataPart->isFCS_gridsize_z()));

	/// Update LineEdit for "Maximum number of iterations". 
	fcsMaxInterationsLineEdit->setText(QString::number(myDataPart->isFCS_max_iterations()));

	/// Update LineEdit for "Debug level". 
	fcsDebuglevelLineEdit->setText(QString::number(myDataPart->isFCS_debuglevel()));

	/// Update LineEdit for "epsilon". 
	fcsEpsilonLineEdit->setText(QString::number(myDataPart->isFCS_epsilon(), 'g'));

	/// Update LineEdit for "NPM". 
	fcsNpmLineEdit->setText(QString::number(myDataPart->isFCS_npm(), 'g'));

	/// Update LineEdit for "Number of walk threads per MPI rank". 
	fcsNumWalkThreadsLineEdit->setText(QString::number(myDataPart->isFCS_num_walk_threads()));

	/// Update LineEdit for "theta". 
	fcsThetaLineEdit->setText(QString::number(myDataPart->isFCS_theta(), 'g'));

	/// Update LineEdit for "Multi-grid cycle type". 
	fcsCycleTypeLineEdit->setText(QString::number(myDataPart->isFCS_cycle_type()));

	/// Update LineEdit for "Order of discretization scheme for PDE solver". 
	fcsDiscretizationOrderLineEdit->setText(QString::number(myDataPart->isFCS_discretization_order()));

	/// Update LineEdit for "Degree of interplation plynomial". 
	fcsInterpolationOrderLineEdit->setText(QString::number(myDataPart->isFCS_interpolation_order()));

	/// Update LineEdit for "Maximum number of multi-grid levels". 
	fcsMaxLevelLineEdit->setText(QString::number(myDataPart->isFCS_max_level()));

	/// Update LineEdit for "Range of smearing of charges over near field cells". 
	fcsNearFieldCellsLineEdit->setText(QString::number(myDataPart->isFCS_near_field_cells()));

	/// Update LineEdit for "Threshold for residual for solver iteration".
	fcsPrecisionLabel->setText(QString::number(myDataPart->isFCS_precision(), 'g'));

	/// Update LineEdit for "Number of smoothing steps". 
	fcsSmoothingStepsLineEdit->setText(QString::number(myDataPart->isFCS_smoothing_steps()));

	/// Update LineEdit for "epsI". 
	fcsEpsILineEdit->setText(QString::number(myDataPart->isFCS_epsI(), 'g'));

	/// Update LineEdit for "m". 
	fcsMLineEdit->setText(QString::number(myDataPart->isFCS_m()));

	/// Update LineEdit for "p". 
	fcsPLineEdit->setText(QString::number(myDataPart->isFCS_p()));

	/// Update LineEdits for "Size of oversampled FFT grid". 
	fcsOversampledGridsizeXLineEdit->setText(QString::number(myDataPart->isFCS_oversampled_gridsize_x()));
	fcsOversampledGridsizeYLineEdit->setText(QString::number(myDataPart->isFCS_oversampled_gridsize_y()));
	fcsOversampledGridsizeZLineEdit->setText(QString::number(myDataPart->isFCS_oversampled_gridsize_z()));

	// Enable SpinBoxes for parallelization cuts. 
	cutsXSpinBox->setEnabled( myDataPart->isSingleProcs() );
	cutsYSpinBox->setEnabled( myDataPart->isSingleProcs() );
	cutsZSpinBox->setEnabled( myDataPart->isSingleProcs() );

	// Update SpinBox for "Used Processes".
	usedProcsSpinBox->setValue( myDataPart->isProcsAtAll() );
	usedProcsValueLabel->setText( QString::number( myDataPart->isProcessorUsage() ) );
	if (myDataPart->isSingleProcs()) usedProcsWidgetStack->raiseWidget( 0 );
	else usedProcsWidgetStack->raiseWidget( 1 );

	// Set values of SpinBoxes for parallelization cuts. 
	cutsXSpinBox->setValue( myDataPart->isCutsX() );
	cutsYSpinBox->setValue( myDataPart->isCutsY() );
	cutsZSpinBox->setValue( myDataPart->isCutsZ() );

	// Updating Widgets for LoadBalancing function. 
	loadBalancingfunctionComboBox->setCurrentText( myDataPart->isLoadBalancingFunction() );
	loadBal_tLineEdit->setText( QString::number( myDataPart->isLoadbal_t(), 'g' ));
	loadBaL_deltaLineEdit->setText( QString::number( myDataPart->isLoadbal_delta(), 'g' ));
	loadBal_stepLineEdit->setText( QString::number( myDataPart->isLoadbal_step() ) );
}

/// this function handles: "linkedCellR_CutLineEdit sends signal: returnpressed!"
void SolParallel_GUI::lCR_CLErp()
{
	myDataPart->setLinkedCellR_Cut( linkedCellR_CutLineEdit->text().toDouble() );
}

/// this function handles: "epsilon0invLineEdit sends signal: returnpressed!"
void SolParallel_GUI::epsilon0inv_rp()
{
	myDataPart->setEpsilon0inv( epsilon0invLineEdit->text().toDouble() );
}

/// this function handles: "r_cutLineEdit sends signal: returnpressed!"
void SolParallel_GUI::R_CLErp()
{
	myDataPart->setR_Cut( r_cutLineEdit->text().toDouble() );
}

/// this function handles: "r_lLineEdit sends signal: returnpressed!"
void SolParallel_GUI::R_LLErp()
{
	myDataPart->setR_L( r_lLineEdit->text().toDouble() );
}

/// this function handles: "splittingCoeffizientGLineEdit sends signal: returnpressed!"
void SolParallel_GUI::sCGLErp()
{
	myDataPart->setSplittingCoefficientG( splittingCoeffizientGLineEdit->text().toDouble() );
}

/// this function handles: "MapLineEdit sends signal: returnpressed!"
void SolParallel_GUI::mapLErp()
{
	myDataPart->setMAP( MapLineEdit->text().toDouble() );
}

/// this function handles: "cellratioLineEdit sends signal: returnpressed!"
void SolParallel_GUI::cellratioLErp()
{
	myDataPart->setCellratio( cellratioLineEdit->text().toInt() );
}

/// this function handles: "interpolationLineEdit sends signal: returnpressed!"
void SolParallel_GUI::interpolationLErp()
{
	myDataPart->setInterpolationDegree( interpolationLineEdit->text().toInt() );
}

/// this function handles: "maxTreeLineEdit sends signal: returnpressed!"
void SolParallel_GUI::maxTreeLErp()
{
	myDataPart->setMaxTreeLevel( maxTreeLineEdit->text().toInt() );
}

/// this function handles "fcsToleranceLineEdit sends singal: returnpressed!"
void SolParallel_GUI::fcstoleranceLErp()
{
	myDataPart->setFCS_tolerance(fcsToleranceLineEdit->text().toDouble());
}

void SolParallel_GUI::fcstoleranceTypeLErp()
{
	myDataPart->setFCS_tolerance_type(fcsToleranceTypeLineEdit->text().toInt());
}

void SolParallel_GUI::fcssCAlphaLErp()
{
	myDataPart->setFCS_splittingCoefficientAlpha(fcsSplittingCoefficientAlphaLineEdit->text().toDouble());
}

void SolParallel_GUI::fcspImagesXLErp()
{
	myDataPart->setFCS_periodic_images_x(fcsPeriodicImagesXLineEdit->text().toInt());
}

void SolParallel_GUI::fcspImagesYLErp()
{
	myDataPart->setFCS_periodic_images_y(fcsPeriodicImagesYLineEdit->text().toInt());
}

void SolParallel_GUI::fcspImagesZLErp()
{
	myDataPart->setFCS_periodic_images_z(fcsPeriodicImagesZLineEdit->text().toInt());
}

void SolParallel_GUI::fcskmaxLErp()
{
	myDataPart->setFCS_kmax(fcsKmaxLineEdit->text().toInt());
}

void SolParallel_GUI::fcsmaxkmaxLErp()
{
	myDataPart->setFCS_maxkmax(fcsMaxkmaxLineEdit->text().toInt());
}

void SolParallel_GUI::fcsbalanceloadLErp()
{
	myDataPart->setFCS_balanceload(fcsBalanceloadLineEdit->text().toInt());
}

void SolParallel_GUI::fcsdipoleCLErp()
{
	myDataPart->setFCS_dipole_correction(fcsDipoleCorrectionLineEdit->text().toInt());
}

void SolParallel_GUI::fcsmaxdepthLErp()
{
	myDataPart->setFCS_maxdepth(fcsMaxdepthLineEdit->text().toLongLong());
}

void SolParallel_GUI::fcspotentialLErp()
{
	myDataPart->setFCS_potential(fcsPotentialLineEdit->text().toInt());
}

void SolParallel_GUI::fcsradiusLErp()
{
	myDataPart->setFCS_radius(fcsRadiusLineEdit->text().toDouble());
}

void SolParallel_GUI::fcsuLimitLErp()
{
	myDataPart->setFCS_unroll_limit(fcsUnrollLimitLineEdit->text().toLongLong());
}

void SolParallel_GUI::fcsdegreeLErp()
{
	myDataPart->setFCS_degree(fcsDegreeLineEdit->text().toInt());
}

void SolParallel_GUI::fcsghostsLErp()
{
	myDataPart->setFCS_ghosts(fcsGhostsLineEdit->text().toInt());
}

void SolParallel_GUI::fcsgSizeXLErp()
{
	myDataPart->setFCS_gridsize_x(fcsGridsizeXLineEdit->text().toInt());
}

void SolParallel_GUI::fcsgSizeYLErp()
{
	myDataPart->setFCS_gridsize_y(fcsGridsizeYLineEdit->text().toInt());
}

void SolParallel_GUI::fcsgSizeZLErp()
{
	myDataPart->setFCS_gridsize_z(fcsGridsizeZLineEdit->text().toInt());
}

void SolParallel_GUI::fcsmIterationsLErp()
{
	myDataPart->setFCS_max_iterations(fcsMaxInterationsLineEdit->text().toInt());
}

void SolParallel_GUI::fcsdebuglevelLErp()
{
	myDataPart->setFCS_debuglevel(fcsDebuglevelLineEdit->text().toInt());
}

void SolParallel_GUI::fcsepsilonLErp()
{
	myDataPart->setFCS_epsilon(fcsEpsilonLineEdit->text().toDouble());
}

void SolParallel_GUI::fcsloadbalancingLErp()
{
	myDataPart->setFCS_load_balancing(fcsLoadBalancingLabel->text().toInt());
}

void SolParallel_GUI::fcsnpmLErp()
{
	myDataPart->setFCS_npm(fcsNpmLineEdit->text().toDouble());
}

void SolParallel_GUI::fcsnwthreadsLErp()
{
	myDataPart->setFCS_num_walk_threads(fcsNumWalkThreadsLineEdit->text().toInt());
}

void SolParallel_GUI::fcsthetaLErp()
{
	myDataPart->setFCS_theta(fcsThetaLineEdit->text().toDouble());
}

void SolParallel_GUI::fcscycletypeLErp()
{
	myDataPart->setFCS_cycle_type(fcsCycleTypeLineEdit->text().toInt());
}

void SolParallel_GUI::fcsdisOrderLErp()
{
	myDataPart->setFCS_discretization_order(fcsDiscretizationOrderLineEdit->text().toInt());
}

void SolParallel_GUI::fcsintOrderLErp()
{
	myDataPart->setFCS_interpolation_order(fcsInterpolationOrderLineEdit->text().toInt());
}

void SolParallel_GUI::fcsmaxLevelLErp()
{
	myDataPart->setFCS_max_level(fcsMaxLevelLineEdit->text().toInt());
}

void SolParallel_GUI::fcsnfcellsLErp()
{
	myDataPart->setFCS_near_field_cells(fcsNearFieldCellsLineEdit->text().toInt());
}

void SolParallel_GUI::fcsprecisionLErp()
{
	myDataPart->setFCS_precision(fcsPrecisionLineEdit->text().toDouble());
}

void SolParallel_GUI::fcssStepsLErp()
{
	myDataPart->setFCS_smoothing_steps(fcsSmoothingStepsLineEdit->text().toInt());
}

void SolParallel_GUI::fcsepsILErp()
{
	myDataPart->setFCS_epsI(fcsEpsILineEdit->text().toDouble());
}

void SolParallel_GUI::fcsmLErp()
{
	myDataPart->setFCS_m(fcsMLineEdit->text().toInt());
}

void SolParallel_GUI::fcspLErp()
{
	myDataPart->setFCS_p(fcsPLineEdit->text().toInt());
}

void SolParallel_GUI::fcsoversampledGsXLErp()
{
	myDataPart->setFCS_oversampled_gridsize_x(fcsOversampledGridsizeXLineEdit->text().toInt());
}

void SolParallel_GUI::fcsoversampledGsYLErp()
{
	myDataPart->setFCS_oversampled_gridsize_y(fcsOversampledGridsizeYLineEdit->text().toInt());
}

void SolParallel_GUI::fcsoversampledGsZLErp()
{
	myDataPart->setFCS_oversampled_gridsize_z(fcsOversampledGridsizeZLineEdit->text().toInt());
}

/// this function handles: "loadBal_tLineEdit sends signal: returnpressed!"
void SolParallel_GUI::lb_tLErp()
{
	myDataPart->setLoadbal_t( loadBal_tLineEdit->text().toDouble() );
}

/// this function handles: "loadBaL_deltaLineEdit sends signal: returnpressed!"
void SolParallel_GUI::lb_deltaLErp()
{
	myDataPart->setLoadbal_delta( loadBaL_deltaLineEdit->text().toDouble() );
}

/// this function handles: "loadBal_stepLineEdit sends signal: returnpressed!"
void SolParallel_GUI::lb_stepLErp()
{
	myDataPart->setLoadbal_step( loadBal_stepLineEdit->text().toInt() );
}

#include "solparallel_gui.moc"
