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

	// Set tolerance widgets. 
	fcsToleranceLabel = new QLabel(DataButtonGroup, "fcsToleranceLabel");
	DataButtonGroupLayout->addWidget(fcsToleranceLabel, 9, 0);
	fcsToleranceLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_tolerance(), 'g'), DataButtonGroup, "fcsToleranceLineEdit");
	QDoubleValidator* fcsToleranceLineEditValidator = new QDoubleValidator(DBL_MIN, DBL_MAX, 8, fcsToleranceLineEdit);
	fcsToleranceLineEdit->setValidator(fcsToleranceLineEditValidator);
	connect(fcsToleranceLineEdit, SIGNAL(returnPressed()), this, SLOT(fcstoleranceLErp()));
	connect(fcsToleranceLineEdit, SIGNAL(lostFocus()), this, SLOT(fcstoleranceLErp()));
	fcsToleranceLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsToleranceLineEdit, 9, 1);
	
	// Set tolerance_type widgets.
	fcsToleranceTypeLabel = new QLabel(DataButtonGroup, "fcsToleranceTypeLabel");
	DataButtonGroupLayout->addWidget(fcsToleranceTypeLabel, 10, 0);
	fcsToleranceTypeLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_tolerance_type()), DataButtonGroup, "fcsToleranceTypeLineEdit");
	QIntValidator* fcsToleranceTypeLineEditValidator = new QIntValidator(0, INT_MAX, fcsToleranceTypeLineEdit);
	fcsToleranceTypeLineEdit->setValidator(fcsToleranceTypeLineEditValidator);
	connect(fcsToleranceTypeLineEdit, SIGNAL(returnPressed()), this, SLOT(fcstoleranceTypeLErp()));
	connect(fcsToleranceTypeLineEdit, SIGNAL(lostFocus()), this, SLOT(fcstoleranceTypeLErp()));
	fcsToleranceTypeLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsToleranceTypeLineEdit, 10, 1);

	// Set periodic_images widgets. 
	fcsPeriodicImagesLabel = new QLabel(DataButtonGroup, "fcsPeriodicImagesLabel");
	DataButtonGroupLayout->addWidget(fcsPeriodicImagesLabel, 11, 0);
	fcsPeriodicImagesXLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_periodic_images_x()), DataButtonGroup, "fcsPeriodicImagesXLineEdit");
	QIntValidator* fcsPeriodicImagesXLineEditValidator = new QIntValidator(0, INT_MAX, fcsPeriodicImagesXLineEdit);
	fcsPeriodicImagesXLineEdit->setValidator(fcsPeriodicImagesXLineEditValidator);
	connect(fcsPeriodicImagesXLineEdit, SIGNAL(returnPressed()), this, SLOT(fcspImagesXLErp()));
	connect(fcsPeriodicImagesXLineEdit, SIGNAL(lostFocus()), this, SLOT(fcspImagesXLErp()));
	fcsPeriodicImagesXLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPeriodicImagesXLineEdit, 11, 1);

	fcsPeriodicImagesYLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_periodic_images_y()), DataButtonGroup, "fcsPeriodicImagesYLineEdit");
	QIntValidator* fcsPeriodicImagesYLineEditValidator = new QIntValidator(0, INT_MAX, fcsPeriodicImagesYLineEdit);
	fcsPeriodicImagesYLineEdit->setValidator(fcsPeriodicImagesYLineEditValidator);
	connect(fcsPeriodicImagesYLineEdit, SIGNAL(returnPressed()), this, SLOT(fcspImagesYLErp()));
	connect(fcsPeriodicImagesYLineEdit, SIGNAL(lostFocus()), this, SLOT(fcspImagesYLErp()));
	fcsPeriodicImagesYLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPeriodicImagesYLineEdit, 11, 2);

	fcsPeriodicImagesZLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_periodic_images_z()), DataButtonGroup, "fcsPeriodicImagesZLineEdit");
	QIntValidator* fcsPeriodicImagesZLineEditValidator = new QIntValidator(0, INT_MAX, fcsPeriodicImagesZLineEdit);
	fcsPeriodicImagesZLineEdit->setValidator(fcsPeriodicImagesZLineEditValidator);
	connect(fcsPeriodicImagesZLineEdit, SIGNAL(returnPressed()), this, SLOT(fcspImagesZLErp()));
	connect(fcsPeriodicImagesZLineEdit, SIGNAL(lostFocus()), this, SLOT(fcspImagesZLErp()));
	fcsPeriodicImagesZLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPeriodicImagesZLineEdit, 11, 3);
	
	// Set kmax widgets. 
	fcsKmaxLabel = new QLabel(DataButtonGroup, "fcsKmaxLabel");
	DataButtonGroupLayout->addWidget(fcsKmaxLabel, 12, 0);
	fcsKmaxLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_kmax()), DataButtonGroup, "fcsKmaxLineEdit");
	QIntValidator* fcsKmaxLineEditValidator = new QIntValidator(0, INT_MAX, fcsKmaxLineEdit);
	fcsKmaxLineEdit->setValidator(fcsKmaxLineEditValidator);
	connect(fcsKmaxLineEdit, SIGNAL(returnPressed()), this, SLOT(fcskmaxLErp()));
	connect(fcsKmaxLineEdit, SIGNAL(lostFocus()), this, SLOT(fcskmaxLErp()));
	fcsKmaxLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsKmaxLineEdit, 12, 1);

	// Set maxkmax widgets. 
	fcsMaxkmaxLabel = new QLabel(DataButtonGroup, "fcsMaxkmaxLabel");
	DataButtonGroupLayout->addWidget(fcsMaxkmaxLabel, 13, 0);
	fcsMaxkmaxLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_maxkmax()), DataButtonGroup, "fcsMaxkmaxLineEdit");
	QIntValidator* fcsMaxkmaxLineEditValidator = new QIntValidator(0, INT_MAX, fcsMaxkmaxLineEdit);
	fcsMaxkmaxLineEdit->setValidator(fcsmaxKmaxLineEditValidator);
	connect(fcsMaxkmaxLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsmaxkmaxLErp()));
	connect(fcsMaxkmaxLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsmaxkmaxLErp()));
	fcsMaxkmaxLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsMaxkmaxLineEdit, 13, 1);

	// Set balanceload widgets.
	fcsBalanceloadLabel = new QLabel(DataButtonGroup, "fcsBalanceloadLabel");
	DataButtonGroupLayout->addWidget(fcsBalanceloadLabel, 14, 0);
	fcsBalanceloadLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_balanceload()), DataButtonGroup, "fcsBalanceloadLineEdit");
	QIntValidator* fcsBalanceloadLineEditValidator = new QIntValidator(0, INT_MAX, fcsBalanceloadLineEdit);
	fcsBalanceloadLineEdit->setValidator(fcsBalanceloadLineEditValidator);
	connect(fcsBalanceloadLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsbalanceloadLErp()));
	connect(fcsBalanceloadLineEdit, SIGNAL(lostFocus()), this, SLOT((fcsbalanceloadLErp())));
	fcsBalanceloadLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsBalanceloadLineEdit, 14, 1);

	// Set dipole_correction widgets.
	fcsDipoleCorrectionLabel = new QLabel(DataButtonGroup, "fcsDipoleCorrectionLabel");
	DataButtonGroupLayout->addWidget(fcsDipoleCorrectionLabel, 15, 0);
	fcsDipoleCorrectionLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_dipole_correction()), DataButtonGroup, "fcsDipoleCorrectionLineEdit");
	QIntValidator* fcsDipoleCorrectionLineEditValidator = new QIntValidator(0, INT_MAX, fcsDipoleCorrectionLineEdit);
	fcsDipoleCorrectionLineEdit->setValidator(fcsDipoleCorrectionLineEditValidator);
	connect(fcsDipoleCorrectionLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsdipoleCLErp()));
	connect(fcsDipoleCorrectionLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsdipoleCLErp()));
	fcsDipoleCorrectionLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsDipoleCorrectionLineEdit, 15, 1);

	// Set maxdepth widgets. 
	fcsMaxdepthLabel= new QLabel(DataButtonGroup, "fcsMaxdepthLabel");
	DataButtonGroupLayout->addWidget(fcsMaxdepthLabel, 16, 0);
	fcsMaxdepthLineEdit= new QLineEdit(QString::number(myDataPart->isFCS_maxdepth()), DataButtonGroup, "fcsMaxdepthLineEdit");
	QIntValidator* fcsMaxdepthLineEditValidator = new QIntValidator(0, INT_MAX, fcsMaxdepthLineEdit);
	fcsMaxdepthLineEdit->setValidator(fcsMaxdepthLineEditValidator);
	connect(fcsMaxdepthLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsmaxdepthLErp()));
	connect(fcsMaxdepthLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsmaxdepthLErp()));
	fcsMaxdepthLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsMaxdepthLineEdit, 16, 1);

	// Set potential widgets. 
	fcsPotentialLabel = new QLabel(DataButtonGroup, "fcsPotentialLabel");
	DataButtonGroupLayout->addWidget(fcsPotentialLabel, 17, 0);
	fcsPotentialLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_potential()), DataButtonGroup, "fcsPotentialLineEdit");
	QIntValidator* fcsPotentialLineEditValidator = new QIntValidator(0, INT_MAX, fcsPotentialLineEdit);
	fcsPotentialLineEdit->setValidator(fcsPotentialLineEditValidator);
	connect(fcsPotentialLineEdit, SIGNAL(returnPressed()), this, SLOT(fcspotentialLErp()));
	connect(fcsPotentialLineEdit, SIGNAL(lostFocus()), this, SLOT(fcspotentialLErp()));
	fcsPotentialLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPotentialLineEdit, 17, 1);

	// Set radius widgets.
	fcsRadiusLabel = new QLabel(DataButtonGroup, "fcsRadiusLabel");
	DataButtonGroupLayout->addWidget(fcsRadiusLabel, 18, 0);
	fcsRadiusLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_radius(), 'g'), DataButtonGroup, "fcsRadiusLineEdit");
	QDoubleValidator* fcsRadiusLineEditValidator = new QDoubleValidator(0, DBL_MAX, 8, fcsRadiusLineEdit);
	fcsRadiusLineEdit->setValidator(fcsRadiusLineEditValidator);
	connect(fcsRadiusLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsradiusLErp()));
	connect(fcsRadiusLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsradiusLErp()));
	fcsRadiusLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsRadiusLineEdit, 18, 1);

	// Set unroll_limit widgets. 
	fcsUnrollLimitLabel = new QLabel(DataButtonGroup, "fcsUnrollLimitLabel");
	DataButtonGroupLayout->addWidget(fcsUnrollLimitLabel, 19, 0);
	fcsUnrollLimitLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_unroll_limit()), DataButtonGroup, "fcsUnrollLimitLineEdit");
	QIntValidator* fcsUnrollLimitLineEditValidator = new QIntValidator(0, INT_MAX, fcsUnrollLimitLineEdit);
	fcsUnrollLimitLineEdit->setValidator(fcsUnrollLimitLineEditValidator);
	connect(fcsUnrollLimitLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsuLimitLErp()));
	connect(fcsUnrollLimitLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsuLimitLErp()));
	fcsUnrollLimitLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsUnrollLimitLineEdit, 19, 1);

	// Set degree widgets. 
	fcsDegreeLabel = new QLabel(DataButtonGroup, "fcsDegreeLabel");
	DataButtonGroupLayout->addWidget(fcsDegreeLabel, 20, 0);
	fcsDegreeLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_degree()), DataButtonGroup, "fcsDegreeLineEdit");
	QIntValidator* fcsDegreeLineEditValidator = new QIntValidator(0, INT_MAX, fcsDegreeLineEdit);
	fcsDegreeLineEdit->setValidator(fcsDegreeLineEditValidator);
	connect(fcsDegreeLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsdegreeLErp()));
	connect(fcsDegreeLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsdegreeLErp()));
	fcsDegreeLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsDegreeLineEdit, 20, 1);

	// Set ghosts widgets. 
	fcsGhostsLabel = new QLabel(DataButtonGroup, "fcsGhostsLabel");
	DataButtonGroupLayout->addWidget(fcsGhostsLabel, 21, 0);
	fcsGhostsLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_ghosts()), DataButtonGroup, "fcsGhostsLineEdit");
	QIntValidator* fcsGhostsLineEditValidator = new QIntValidator(0, INT_MAX, fcsGhostsLineEdit);
	fcsGhostsLineEdit->setValidator(fcsGhostsLineEditValidator);
	connect(fcsGhostsLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsghostsLErp()));
	connect(fcsGhostsLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsghostsLErp()));
	fcsGhostsLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsGhostsLineEdit, 21, 1);

	// Set gridsize widgets. 
	fcsGridsizeLabel = new QLabel(DataButtonGroup, "fcsGridsizeLabel");
	DataButtonGroupLayout->addWidget(fcsGridsizeLabel, 22, 0);
	// x-axis
	fcsGridsizeXLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_gridsize_x()), DataButtonGroup, "fcsGridsizeXLineEdit");
	QIntValidator* fcsGridsizeXLineEditValidator = new QIntValidator(0, INT_MAX, fcsGridsizeXLineEdit);
	fcsGridsizeXLineEdit->setValidator(fcsGridsizeXLineEditValidator);
	connect(fcsGridsizeXLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsgSizeXLErp()));
	connect(fcsGridsizeXLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsgSizeXLErp()));
	fcsGridsizeXLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsGridsizeXLineEdit, 22, 1);
	// y-axis
	fcsGridsizeYLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_gridsize_y()), DataButtonGroup, "fcsGridsizeYLineEdit");
	QIntValidator* fcsGridsizeYLineEditValidator = new QIntValidator(0, INT_MAX, fcsGridsizeYLineEdit);
	fcsGridsizeYLineEdit->setValidator(fcsGridsizeYLineEditValidator);
	connect(fcsGridsizeYLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsgSizeYLErp()));
	connect(fcsGridsizeYLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsgSizeYLErp()));
	fcsGridsizeYLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsGridsizeYLineEdit, 22, 2);
	// z-axis
	fcsGridsizeZLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_gridsize_z()), DataButtonGroup, "fcsGridsizeZLineEdit");
	QIntValidator* fcsGridsizeZLineEditValidator = new QIntValidator(0, INT_MAX, fcsGridsizeZLineEdit);
	fcsGridsizeZLineEdit->setValidator(fcsGridsizeZLineEditValidator);
	connect(fcsGridsizeZLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsgSizeZLErp()));
	connect(fcsGridsizeZLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsgSizeZLErp()));
	fcsGridsizeZLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsGridsizeZLineEdit, 22, 3);

	// Set max_iterations widgets. 
	fcsMaxIterationsLabel = new QLabel(DataButtonGroup, "fcsMaxIterationsLabel");
	DataButtonGroupLayout->addWidget(fcsMaxIterationsLabel, 23, 0);
	fcsMaxIterationsLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_max_iterations()), DataButtonGroup, "fcsMaxIterationsLineEdit");
	QIntValidator* fcsMaxIterationsLineEditValidator = new QIntValidator(0, INT_MAX, fcsMaxIterationsLineEdit);
	fcsMaxIterationsLineEdit->setValidator(fcsMaxIterationsLineEditValidator);
	connect(fcsMaxIterationsLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsmIterationsLErp()));
	connect(fcsMaxIterationsLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsmIterationsLErp()));
	fcsMaxIterationsLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsMaxIterationsLineEdit, 23, 1);

	// Set debuglevel widgets. 
	fcsDebuglevelLabel = new QLabel(DataButtonGroup, "fcsDebuglevelLabel");
	DataButtonGroupLayout->addWidget(fcsDebuglevelLabel, 24, 0);
	fcsDebuglevelLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_debuglevel()), DataButtonGroup, "fcsDebuglevelLineEdit");
	QIntValidator* fcsDebuglevelLineEditValidator = new QIntValidator(0, INT_MAX, fcsDebuglevelLineEdit);
	fcsDebuglevelLineEdit->setValidator(fcsDebuglevelLineEditValidator);
	connect(fcsDebuglevelLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsdebuglevelLErp()));
	connect(fcsDebuglevelLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsdebuglevelLErp()));
	fcsDebuglevelLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsDebuglevelLineEdit, 24, 1);

	// Set epsilon widgets. 
	fcsEpsilonLabel = new QLabel(DataButtonGroup, "fcsEpsilonLabel");
	DataButtonGroupLayout->addWidget(fcsEpsilonLabel, 25, 0);
	fcsEpsilonLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_epsilon(), 'g'), DataButtonGroup, "fcsEpsilonLineEdit");
	QDoubleValidator* fcsEpsilonLineEditValidator = new QDoubleValidator(0, DBL_MAX, 8, fcsEpsilonLineEdit);
	fcsEpsilonLineEdit->setValidator(fcsEpsilonLineEditValidator);
	connect(fcsEpsilonLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsepsilonLErp()));
	connect(fcsEpsilonLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsepsilonLErp()));
	fcsEpsilonLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsEpsilonLineEdit, 25, 1);

	// Set load_balancing widgets.
	fcsLoadBalancingLabel = new QLabel(DataButtonGroup, "fcsLoadBalancingLabel");
	DataButtonGroupLayout->addWidget(fcsLoadBalancingLabel, 26, 0);
	fcsLoadBalancingLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_load_balancing()), DataButtonGroup, "fcsLoadBalancingLineEdit");
	QIntValidator* fcsLoadBalancingLineEditValidator = new QIntValidator(0, INT_MAX, fcsLoadBalancingLineEdit);
	fcsLoadBalancingLineEdit->setValidator(fcsLoadBalancingLineEditValidator);
	connect(fcsLoadBalancingLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsloadbalancingLErp()));
	connect(fcsLoadBalancingLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsloadbalancingLErp()));
	fcsLoadBalancingLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsLoadBalancingLineEdit, 26, 1);

	// Set npm widgets. 
	fcsNpmLabel = new QLabel(DataButtonGroup, "fcsNpmLabel");
	DataButtonGroupLayout->addWidget(fcsNpmLabel, 27, 0);
	fcsNpmLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_npm(), 'g'), DataButtonGroup, "fcsNpmLineEdit");
	QDoubleValidator* fcsNpmLineEditValidator = new QDoubleValidator(DBL_MIN, DBL_MAX, 8, fcsNpmLineEdit);
	fcsNpmLineEdit->setValidator(fcsNpmLineEditValidator);
	connect(fcsNpmLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsnpmLErp()));
	connect(fcsNpmLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsnpmLErp()));
	fcsNpmLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsNpmLineEdit, 27, 1);

	// Set num_walk_threads widgets. 
	fcsNumWalkThreadsLabel = new QLabel(DataButtonGroup, "fcsNumWalkThreadsLabel");
	DataButtonGroupLayout->addWidget(fcsNumWalkThreadsLabel, 28, 0);
	fcsNumWalkThreadsLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_num_walk_threads()), DataButtonGroup, "fcsNumWalkThreadsLineEdit");
	QIntValidator* fcsNumWalkThreadsLineEditValidator = new QIntValidator(0, INT_MAX, fcsNumWalkThreadsLineEdit);
	fcsNumWalkThreadsLineEdit->setValidator(fcsNumWalkThreadsLineEditValidator);
	connect(fcsNumWalkThreadsLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsnwthreadsLErp()));
	connect(fcsNumWalkThreadsLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsnwthreadsLErp()));
	fcsNumWalkThreadsLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsNumWalkThreadsLineEdit, 28, 1);

	// Set theta widgets. 
	fcsThetaLabel = new QLabel(DataButtonGroup, "fcsThetaLabel");
	DataButtonGroupLayout->addWidget(fcsThetaLabel, 29, 0);
	fcsThetaLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_theta(), 'g'), DataButtonGroup, "fcsThetaLineEdit");
	QDoubleValidator* fcsThetaLineEditValidator = new QDoubleValidator(DBL_MIN, DBL_MAX, 8, fcsThetaLineEdit);
	fcsThetaLineEdit->setValidator(fcsThetaLineEditValidator);
	connect(fcsThetaLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsthetaLErp()));
	connect(fcsThetaLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsthetaLErp()));
	fcsThetaLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsThetaLineEdit, 29, 1);

	// Set cycle_type widgets. 
	fcsCycleTypeLabel = new QLabel(DataButtonGroup, "fcsCycleTypeLabel");
	DataButtonGroupLayout->addWidget(fcsCycleTypeLabel, 30, 0);
	fcsCycleTypeLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_cycle_type()), DataButtonGroup, "fcsCycleTypeLineEdit");
	QIntValidator* fcsCycleTypeLineEditValidator = new QIntValidator(0, INT_MAX, fcsCycleTypeLineEdit);
	fcsCycleTypeLineEdit->setValidator(fcsCycleTypeLineEditValidator);
	connect(fcsCycleTypeLineEdit, SIGNAL(returnPressed()), this, SLOT(fcscycletypeLErp()));
	connect(fcsCycleTypeLineEdit, SIGNAL(lostFocus()), this, SLOT(fcscycletypeLErp()));
	fcsCycleTypeLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsCycleTypeLineEdit, 30, 1);

	// Set discretization_order widgets. 
	fcsDiscretizationOrderLabel = new QLabel(DataButtonGroup, "fcsDiscretizationOrderLabel");
	DataButtonGroupLayout->addWidget(fcsDiscretizationOrderLabel, 31, 0);
	fcsDiscretizationOrderLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_discretization_order()), DataButtonGroup, "fcsDiscretizationOrderLineEdit");
	QIntValidator* fcsDiscretizationOrderLineEditValidator = new QIntValidator(0, INT_MAX, fcsDiscretizationOrderLineEdit);
	fcsDiscretizationOrderLineEdit->setValidator(fcsDiscretizationOrderLineEditValidator);
	connect(fcsDiscretizationOrderLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsdisOrderLErp()));
	connect(fcsDiscretizationOrderLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsdisOrderLErp()));
	fcsDiscretizationOrderLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsDiscretizationOrderLineEdit, 31, 1);

	// Set interpolation_order widgets. 
	fcsInterpolationOrderLabel = new QLabel(DataButtonGroup, "fcsInterpolationOrderLabel");
	DataButtonGroupLayout->addWidget(fcsInterpolationOrderLabel, 32, 0);
	fcsInterpolationOrderLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_interpolation_order()), DataButtonGroup, "fcsInterpolationOrderLineEdit");
	QIntValidator* fcsInterpolationOrderLineEditValidator = new QIntValidator(0, INT_MAX, fcsInterpolationOrderLineEdit);
	fcsInterpolationOrderLineEdit->setValidator(fcsInterpolationOrderLineEditValidator);
	connect(fcsInterpolationOrderLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsintOrderLErp()));
	connect(fcsInterpolationOrderLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsintOrderLErp()));
	fcsInterpolationOrderLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsInterpolationOrderLineEdit, 32, 1);

	// Set max_level widgets. 
	fcsMaxLevelLabel = new QLabel(DataButtonGroup, "fcsMaxLevelLabel");
	DataButtonGroupLayout->addWidget(fcsMaxLevelLabel, 33, 0);
	fcsMaxLevelLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_max_level()), DataButtonGroup, "fcsMaxLevelLineEdit");
	QIntValidator* fcsMaxLevelLineEditValidator = new QIntValidator(0, INT_MAX, fcsMaxLevelLineEdit);
	fcsMaxLevelLineEdit->setValidator(fcsMaxLevelLineEditValidator);
	connect(fcsMaxLevelLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsMaxLevelLineEdit));
	connect(fcsMaxLevelLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsMaxLevelLineEdit));
	fcsMaxLevelLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsMaxLevelLineEdit, 33, 1);

	// Set near_field_cells widgets. 
	fcsNearFieldCellsLabel = new QLabel(DataButtonGroup, "fcsNearFieldCellsLabel");
	DataButtonGroupLayout->addWidget(fcsNearFieldCellsLabel, 34, 0);
	fcsNearFieldCellsLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_near_field_cells()), DataButtonGroup, "fcsNearFieldCellsLineEdit");
	QIntValidator* fcsNearFieldCellsLineEditValidator = new QIntValidator(0, INT_MAX, fcsNearFieldCellsLineEdit);
	fcsNearFieldCellsLineEdit->setValidator(fcsNearFieldCellsLineEditValidator);
	connect(fcsNearFieldCellsLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsnfcellsLErp()));
	connect(fcsNearFieldCellsLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsnfcellsLErp()));
	fcsNearFieldCellsLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsNearFieldCellsLineEdit, 34, 1);

	// Set precision widgets. 
	fcsPrecisionLabel = new QLabel(DataButtonGroup, "fcsPrecisionLabel");
	DataButtonGroupLayout->addWidget(fcsPrecisionLabel, 35, 0);
	fcsPrecisionLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_precision(), 'g'), DataButtonGroup, "fcsPrecisionLineEdit");
	QDoubleValidator* fcsPrecisionLineEditValidator = new QDoubleValidator(0, DBL_MAX, fcsPrecisionLineEdit);
	fcsPrecisionLineEdit->setValidator(fcsPrecisionLineEditValidator);
	connect(fcsPrecisionLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsprecisionLErp()));
	connect(fcsPrecisionLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsprecisionLErp()));
	fcsPrecisionLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPrecisionLineEdit, 35, 1);

	// Set smoothing_steps widgets. 
	fcsSmoothingStepsLabel = new QLabel(DataButtonGroup, "fcsSmoothingStepsLabel");
	DataButtonGroupLayout->addWidget(fcsSmoothingStepsLabel, 36, 0);
	fcsSmoothingStepsLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_smoothing_steps()), DataButtonGroup, "fcsSmoothingStepsLineEdit");
	QIntValidator* fcsSmoothingStepsLineEditValidator = new QIntValidator(0, INT_MAX, fcsSmoothingStepsLineEdit);
	fcsSmoothingStepsLineEdit->setValidator(fcsSmoothingStepsLineEditValidator);
	connect(fcsSmoothingStepsLineEdit, SIGNAL(returnPressed()), this, SLOT(fcssStepsLErp()));
	connect(fcsSmoothingStepsLineEdit, SIGNAL(lostFocus()), this, SLOT(fcssStepsLErp()));
	fcsSmoothingStepsLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsSmoothingStepsLineEdit, 36, 1);

	// Set epsI widgets. 
	fcsEpsILabel= new QLabel(DataButtonGroup, "fcsEpsILabel");
	DataButtonGroupLayout->addWidget(fcsEpsILabel, 37, 0);
	fcsEpsILineEdit= new QLineEdit(QString::number(myDataPart->isFCS_epsI(), 'g'), DataButtonGroup, "fcsEpsILineEdit");
	QDoubleValidator* fcsEpsILineEditValidator = new QDoubleValidator(0, DBL_MAX, 8, fcsEpsILineEdit);
	fcsEpsILineEdit->setValidator(fcsEpsILineEditValidator);
	connect(fcsEpsILineEdit, SIGNAL(returnPressed()), this, SLOT(fcsepsILErp()));
	connect(fcsEpsILineEdit, SIGNAL(lostFocus()), this, SLOT(fcsepsILErp()));
	fcsEpsILineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsEpsILineEdit, 37, 1);

	// Set m widgets. 
	fcsMLabel = new QLabel(DataButtonGroup, "fcsMLabel");
	DataButtonGroupLayout->addWidget(fcsMLabel, 38, 0);
	fcsMLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_m()), DataButtonGroup, "fcsMLineEdit");
	QIntValidator* fcsMLineEditValidator = new QIntValidator(INT_MIN, INT_MAX, fcsMLineEdit);
	fcsMLineEdit->setValidator(fcsMLineEditValidator);
	connect(fcsMLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsmLErp()));
	connect(fcsMLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsmLErp()));
	fcsMLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsMLineEdit, 38, 1);

	// Set p widgets. 
	fcsPLabel = new QLabel(DataButtonGroup, "fcsPLabel");
	DataButtonGroupLayout->addWidget(fcsPLabel, 39, 0);
	fcsPLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_p()), DataButtonGroup, "fcsPLineEdit");
	QIntValidator* fcsPLineEditValidator = new QIntValidator(INT_MIN, INT_MAX, fcsPLineEdit);
	fcsPLineEdit->setValidator(fcsPLineEditValidator);
	connect(fcsPLineEdit, SIGNAL(returnPressed()), this, SLOT(fcspLErp()));
	connect(fcsPLineEdit, SIGNAL(lostFocus()), this, SLOT(fcspLErp()));
	fcsPLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsPLineEdit, 39, 1);

	// Set oversampled_gridsize widgets. 
	fcsOversampledGridsizeLabel = new QLabel(DataButtonGroup, "fcsOversampledGridsizeLabel");
	DataButtonGroupLayout->addWidget(fcsOversampledGridsizeLabel, 40, 0);
	// x-axis
	fcsOversampledGridsizeXLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_oversampled_gridsize_x()), DataButtonGroup, "fcsOversampledGridsizeXLineEdit");
	QIntValidator* fcsOversampledGridsizeXLineEditValidator = new QIntValidator(0, INT_MAX, "fcsOversampledGridsizeXLineEdit");
	fcsOversampledGridsizeXLineEdit->setValidator(fcsOversampledGridsizeXLineEditValidator);
	connect(fcsOversampledGridsizeXLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsoversampledGsXLErp()));
	connect(fcsOversampledGridsizeXLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsoversampledGsXLErp()));
	fcsOversampledGridsizeXLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsOversampledGridsizeXLineEdit, 40, 1);
	// y-axis
	fcsOversampledGridsizeYLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_oversampled_gridsize_y()), DataButtonGroup, "fcsOversampledGridsizeYLineEdit");
	QIntValidator* fcsOversampledGridsizeYLineEditValidator = new QIntValidator(0, INT_MAX, "fcsOversampledGridsizeYLineEdit");
	fcsOversampledGridsizeYLineEdit->setValidator(fcsOversampledGridsizeYLineEditValidator);
	connect(fcsOversampledGridsizeYLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsoversampledGsYLErp()));
	connect(fcsOversampledGridsizeYLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsoversampledGsYLErp()));
	fcsOversampledGridsizeYLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsOversampledGridsizeYLineEdit, 40, 1);
	// z-axis
	fcsOversampledGridsizeZLineEdit = new QLineEdit(QString::number(myDataPart->isFCS_oversampled_gridsize_z()), DataButtonGroup, "fcsOversampledGridsizeZLineEdit");
	QIntValidator* fcsOversampledGridsizeZLineEditValidator = new QIntValidator(0, INT_MAX, "fcsOversampledGridsizeZLineEdit");
	fcsOversampledGridsizeZLineEdit->setValidator(fcsOversampledGridsizeZLineEditValidator);
	connect(fcsOversampledGridsizeZLineEdit, SIGNAL(returnPressed()), this, SLOT(fcsoversampledGsZLErp()));
	connect(fcsOversampledGridsizeZLineEdit, SIGNAL(lostFocus()), this, SLOT(fcsoversampledGsZLErp()));
	fcsOversampledGridsizeZLineEdit->setEnabled(false);
	DataButtonGroupLayout->addWidget(fcsOversampledGridsizeZLineEdit, 40, 1);

     
	// Configure visibility and enabling of poisson solver widgets.
	connect(myDataPart, SIGNAL( enablepoisson_solver( bool ) ), poisson_solverComboBox, SLOT( setEnabled( bool ) ) );
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

	// Configure visibility and enabling of cellratio widgets. 
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLineEdit, SLOT( setEnabled( bool ) ) );
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLineEdit, SLOT( setVisible( bool ) ) );
	connect( myDataPart, SIGNAL( enableCellratio( bool ) ), cellratioLabel, SLOT( setVisible( bool ) ) );

	// Configure visibility and enabling of tolerance widgets.
	connect(myDataPart, SIGNAL(enableFCS_tolerance(bool)), fcsToleranceLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_tolerance(bool)), fcsToleranceLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_tolerance(bool)), fcsToleranceLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of tolerance_type widgets. 
	connect(myDataPart, SIGNAL(enableFCS_tolerance_type(bool)), fcsToleranceTypeLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_tolerance_type(bool)), fcsToleranceTypeLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_tolerance_type(bool)), fcsToleranceTypeLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of periodic_images widgets. 
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_x(bool)), fcsPeriodicImagesXLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_x(bool)), fcsPeriodicImagesXLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_y(bool)), fcsPeriodicImagesYLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_y(bool)), fcsPeriodicImagesYLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_z(bool)), fcsPeriodicImagesZLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_z(bool)), fcsPeriodicImagesZLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_x(bool)), fcsPeriodicImagesLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_y(bool)), fcsPeriodicImagesLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_periodic_images_z(bool)), fcsPeriodicImagesLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of splittingCoefficientAlpha widgets.
	connect(myDataPart, SIGNAL(enableFCS_splittingCoefficientAlpha(bool)), fcsSplittingCoefficientAlphaLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_splittingCoefficientAlpha(bool)), fcsSplittingCoefficientAlphaLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_splittingCoefficientAlpha(bool)), fcsSplittingCoefficientAlphaLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of kmax widgets. 
	connect(myDataPart, SIGNAL(enableFCS_kmax(bool)), fcsKmaxLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_kmax(bool)), fcsKmaxLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_kmax(bool)), fcsKmaxLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of maxkmax widgets. 
	connect(myDataPart, SIGNAL(enableFCS_maxkmax(bool)), fcsMaxkmaxLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_maxkmax(bool)), fcsMaxkmaxLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_maxkmax(bool)), fcsMaxkmaxLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of balanceload widgets. 
	connect(myDataPart, SIGNAL(enableFCS_balanceload(bool)), fcsBalanceloadLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_balanceload(bool)), fcsBalanceloadLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_balanceload(bool)), fcsBalanceloadLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of dipole_correction widgets. 
	connect(myDataPart, SIGNAL(enableFCS_dipole_correction(bool)), fcsDipoleCorrectionLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_dipole_correction(bool)), fcsDipoleCorrectionLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_dipole_correction(bool)), fcsDipoleCorrectionLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of maxdepts widgets. 
	connect(myDataPart, SIGNAL(enableFCS_maxdepth(bool)), fcsMaxdepthLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_maxdepth(bool)), fcsMaxdepthLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_maxdepth(bool)), fcsMaxdepthLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of potential widgets. 
	connect(myDataPart, SIGNAL(enableFCS_potential(bool)), fcsPotentialLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_potential(bool)), fcsPotentialLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_potential(bool)), fcsPotentialLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of radius widgets. 
	connect(myDataPart, SIGNAL(enableFCS_radius(bool)), fcsRadiusLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_radius(bool)), fcsRadiusLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_radius(bool)), fcsRadiusLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of unroll_limit widgets. 
	connect(myDataPart, SIGNAL(enableFCS_unroll_limit(bool)), fcsUnrollLimitLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_unroll_limit(bool)), fcsUnrollLimitLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_unroll_limit(bool)), fcsUnrollLimitLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of degree widgets. 
	connect(myDataPart, SIGNAL(enableFCS_degree(bool)), fcsDegreeLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_degree(bool)), fcsDegreeLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_degree(bool)), fcsDegreeLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of ghosts widgets. 
	connect(myDataPart, SIGNAL(enableFCS_ghosts(bool)), fcsGhostsLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_ghosts(bool)), fcsGhostsLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_ghosts(bool)), fcsGhostsLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of gridsize widgets. 
	connect(myDataPart, SIGNAL(enableFCS_gridsize_x(bool)), fcsGridsizeXLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_x(bool)), fcsGridsizeXLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_x(bool)), fcsGridsizeLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_y(bool)), fcsGridsizeYLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_y(bool)), fcsGridsizeYLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_y(bool)), fcsGridsizeLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_z(bool)), fcsGridsizeZLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_z(bool)), fcsGridsizeZLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_gridsize_z(bool)), fcsGridsizeLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of max_iterations widgets. 
	connect(myDataPart, SIGNAL(enableFCS_max_iterations(bool)), fcsMaxIterationsLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_max_iterations(bool)), fcsMaxIterationsLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_max_iterations(bool)), fcsMaxIterationsLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of debuglevel widgets. 
	connect(myDataPart, SIGNAL(enableFCS_debuglevel(bool)), fcsDebuglevelLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_debuglevel(bool)), fcsDebuglevelLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_debuglevel(bool)), fcsDebuglevelLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of epsilon widgets. 
	connect(myDataPart, SIGNAL(enableFCS_epsilon(bool)), fcsEpsilonLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_epsilon(bool)), fcsEpsilonLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_epsilon(bool)), fcsEpsilonLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of load_balancing widgets. 
	connect(myDataPart, SIGNAL(enableFCS_load_balancing(bool)), fcsLoadBalancingLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_load_balancing(bool)), fcsLoadBalancingLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_load_balancing(bool)), fcsLoadBalancingLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of NPM widgets. 
	connect(myDataPart, SIGNAL(enableFCS_npm(bool)), fcsNpmLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_npm(bool)), fcsNpmLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_npm(bool)), fcsNpmLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of num_walk_threads widgets. 
	connect(myDataPart, SIGNAL(enableFCS_num_walk_threads(bool)), fcsNumWalkThreadsLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_num_walk_threads(bool)), fcsNumWalkThreadsLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_num_walk_threads(bool)), fcsNumWalkThreadsLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of theta widgets. 
	connect(myDataPart, SIGNAL(enableFCS_theta(bool)), fcsThetaLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_theta(bool)), fcsThetaLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_theta(bool)), fcsThetaLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of cycle_type widgets. 
	connect(myDataPart, SIGNAL(enableFCS_cycle_type(bool)), fcsCycleTypeLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_cycle_type(bool)), fcsCycleTypeLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_cycle_type(bool)), fcsCycleTypeLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of discretization_order widgets. 
	connect(myDataPart, SIGNAL(enableFCS_discretization_order(bool)), fcsDiscretizationOrderLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_discretization_order(bool)), fcsDiscretizationOrderLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_discretization_order(bool)), fcsDipoleCorrectionLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of interpolation_order widgers. 
	connect(myDataPart, SIGNAL(enableFCS_interpolation_order(bool)), fcsInterpolationOrderLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_interpolation_order(bool)), fcsInterpolationOrderLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_interpolation_order(bool)), fcsInterpolationOrderLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of max_level widgets. 
	connect(myDataPart, SIGNAL(enableFCS_max_level(bool)), fcsMaxLevelLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_max_level(bool)), fcsMaxLevelLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_max_level(bool)), fcsMaxLevelLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of near_field_cells widgets. 
	connect(myDataPart, SIGNAL(enableFCS_near_field_cells(bool)), fcsNearFieldCellsLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_near_field_cells(bool)), fcsNearFieldCellsLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_near_field_cells(bool)), fcsNearFieldCellsLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of precision widgets. 
	connect(myDataPart, SIGNAL(enableFCS_precision(bool)), fcsPrecisionLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_precision(bool)), fcsPrecisionLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_precision(bool)), fcsPrecisionLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of smoothing_steps widgets. 
	connect(myDataPart, SIGNAL(enableFCS_smoothing_steps(bool)), fcsSmoothingStepsLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_smoothing_steps(bool)), fcsSmoothingStepsLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_smoothing_steps(bool)), fcsSmoothingStepsLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of epsI widgets. 
	connect(myDataPart, SIGNAL(enableFCS_epsI(bool)), fcsEpsILineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_epsI(bool)), fcsEpsILineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_epsI(bool)), fcsEpsILabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of m widgets. 
	connect(myDataPart, SIGNAL(enableFCS_m(bool)), fcsMLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_m(bool)), fcsMLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_m(bool)), fcsMLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of p widgets.
	connect(myDataPart, SIGNAL(enableFCS_p(bool)), fcsPLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_p(bool)), fcsPLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_p(bool)), fcsPLabel, SLOT(setVisible(bool)));

	// Configure visibility and enabling of oversampled_gridsize widgets. 
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_x(bool)), fcsOversampledGridsizeXLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_x(bool)), fcsOversampledGridsizeXLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_x(bool)), fcsOversampledGridsizeLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_y(bool)), fcsOversampledGridsizeYLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_y(bool)), fcsOversampledGridsizeYLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_y(bool)), fcsOversampledGridsizeLabel, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_z(bool)), fcsOversampledGridsizeZLineEdit, SLOT(setEnabled(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_z(bool)), fcsOversampledGridsizeZLineEdit, SLOT(setVisible(bool)));
	connect(myDataPart, SIGNAL(enableFCS_oversampled_gridsize_z(bool)), fcsOversampledGridsizeLabel, SLOT(setVisible(bool)));


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
