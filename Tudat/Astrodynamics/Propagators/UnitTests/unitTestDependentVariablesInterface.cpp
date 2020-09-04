/* git    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/Astrodynamics/Propagators/dependentVariablesInterface.h"


namespace tudat
{

namespace unit_tests
{

//Using declarations.
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_dependent_variables_interface )

//template< typename TimeType = double , typename StateScalarType  = double >
//        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
//executeEarthMoonSimulation(
//        const std::vector< std::string > centralBodies,
//        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference =
//        Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
//        const int propagationType = 0,
//        const Eigen::Vector3d parameterPerturbation = Eigen::Vector3d::Zero( ),
//        const bool propagateVariationalEquations = 1 )
//{

//    //Load spice kernels.
//    spice_interface::loadStandardSpiceKernels( );

//    // Define
//    std::vector< std::string > bodyNames;
//    bodyNames.push_back( "Earth" );
//    bodyNames.push_back( "Sun" );
//    bodyNames.push_back( "Moon" );
//    bodyNames.push_back( "Mars" );

//    // Specify initial time
//    TimeType initialEphemerisTime = TimeType( 1.0E7 );
//    TimeType finalEphemerisTime = initialEphemerisTime + 0.5E7;
//    double maximumTimeStep = 3600.0;

//    double buffer = 10.0 * maximumTimeStep;

//    // Create bodies needed in simulation
//    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
//            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

//    NamedBodyMap bodyMap = createBodies( bodySettings );
//    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


//    // Set accelerations between bodies that are to be taken into account.
//    SelectedAccelerationMap accelerationMap;
//    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
//    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    accelerationMap[ "Earth" ] = accelerationsOfEarth;

//    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
//    accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    accelerationMap[ "Moon" ] = accelerationsOfMoon;

//    // Set bodies for which initial state is to be estimated and integrated.
//    std::vector< std::string > bodiesToIntegrate;
//    bodiesToIntegrate.push_back( "Moon" );
//    bodiesToIntegrate.push_back( "Earth" );

//    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

//    // Define propagator settings.
//    std::map< std::string, std::string > centralBodyMap;

//    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
//    {
//        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
//    }

//    // Create acceleration models
//    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//                bodyMap, accelerationMap, centralBodyMap );

//    // Create integrator settings
//    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
//            std::make_shared< IntegratorSettings< TimeType > >
//            ( rungeKutta4, TimeType( initialEphemerisTime ), 1800.0 );


//    // Set initial states of bodies to integrate.
//    TimeType initialIntegrationTime = initialEphemerisTime;

//    // Set (perturbed) initial state.
//    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialTranslationalState;
//    initialTranslationalState = getInitialStatesOfBodies< TimeType, StateScalarType >(
//                bodiesToIntegrate, centralBodies, bodyMap, initialIntegrationTime );
//    initialTranslationalState += initialStateDifference;

//    // Create propagator settings
//    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings;
//    TranslationalPropagatorType propagatorType;
//    if( propagationType == 0 )
//    {
//        propagatorType = cowell;
//    }
//    else if( propagationType == 1 )
//    {
//        propagatorType = encke;
//    }
//    propagatorSettings =  std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
//            ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialTranslationalState,
//              TimeType( finalEphemerisTime ), propagatorType );

//    // Define parameters.
//    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//    {
//        parameterNames.push_back(
//                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
//                        "Moon", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
//                            "Moon", centralBodies[ 0 ], bodyMap, TimeType( initialEphemerisTime ) ) +
//                    initialStateDifference.segment( 0, 6 ),
//                    centralBodies[ 0 ] ) );
//        parameterNames.push_back(
//                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
//                        "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
//                            "Earth", centralBodies[ 1 ], bodyMap, TimeType( initialEphemerisTime ) ) +
//                    initialStateDifference.segment( 6, 6 ),
//                    centralBodies[ 1 ] ) );
//        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );
//        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter ) );
//        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );

//    }

//    // Create parameters
//    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
//            createParametersToEstimate( parameterNames, bodyMap );

//    // Perturb parameters.
//    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
//            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
//    parameterVector.block( 12, 0, 3, 1 ) += parameterPerturbation;
//    parametersToEstimate->resetParameterValues( parameterVector );

//    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

//    {
//        // Create dynamics simulator
//        SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
//                SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
//                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
//                    1, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 1, 0 );

//        // Propagate requested equations.
//        if( propagateVariationalEquations )
//        {
//            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
//        }
//        else
//        {
//            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
//        }

//        // Retrieve test data
//        double testEpoch = 1.4E7;
//        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates =
//                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
//        testStates.block( 0, 0, 6, 1 ) = bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( testEpoch );
//        testStates.block( 6, 0, 6, 1 ) = bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( testEpoch );

//        if( propagateVariationalEquations )
//        {
//            results.first.push_back( dynamicsSimulator.getStateTransitionMatrixInterface( )->
//                                     getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
//        }
//        results.second.push_back( testStates );
//    }
//    return results;
//}

//! Test the dependent variables computation against their numerical propagation.
/*!
 *  Test the dependent variables computation against their numerical propagation. This unit test
 *  ...
 */
BOOST_AUTO_TEST_CASE( testSingleArcDependentVariablesInterface )
{

//    // Load Spice kernels.
//    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;


    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );


    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "ECLIPJ2000" );
    }
    bodySettings[ "Phobos" ] = std::make_shared< BodySettings >( );
    bodySettings[ "Phobos" ]->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
    bodySettings[ "Phobos" ]->gravityFieldSettings =
            getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    NamedBodyMap bodyMap = createBodies( bodySettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "AlienSpaceship" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "AlienSpaceship" ]->setConstantBodyMass( 400.0 );

    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), "Phobos", "ECLIPJ2000" ) );

//    // Create aerodynamic coefficient interface settings.
//    double referenceArea = 4.0;
//    double aerodynamicCoefficient = 1.2;
//    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
//            std::make_shared< ConstantAerodynamicCoefficientSettings >(
//                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

//    // Create and set aerodynamic coefficients object
//    bodyMap[ "AlienSpaceship" ]->setAerodynamicCoefficientInterface(
//                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "AlienSpaceship" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > alienSpaceshipRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "AlienSpaceship" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    alienSpaceshipRadiationPressureSettings, "AlienSpaceship", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAlienSpaceship;
    accelerationsOfAlienSpaceship[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAlienSpaceship[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAlienSpaceship[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationMap[ "AlienSpaceship" ] = accelerationsOfAlienSpaceship;
    bodiesToPropagate.push_back( "AlienSpaceship" );
    centralBodies.push_back( "Phobos" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d alienSpaceshipInitialStateInKeplerianElements;
    alienSpaceshipInitialStateInKeplerianElements( semiMajorAxisIndex ) = 100.0E3;
    alienSpaceshipInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
    alienSpaceshipInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
    alienSpaceshipInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    alienSpaceshipInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    alienSpaceshipInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    double phobosGravitationalParameter = bodyMap.at( "Phobos" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d alienSpaceshipInitialState = convertKeplerianToCartesianElements(
                alienSpaceshipInitialStateInKeplerianElements, phobosGravitationalParameter );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;


    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_acceleration_dependent_variable, "AlienSpaceship" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Sun", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Mars", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Phobos", point_mass_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Earth", third_body_spherical_harmonic_gravity, "AlienSpaceship", "Phobos" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
//                    "AlienSpaceship", "Sun", cannon_ball_radiation_pressure, "AlienSpaceship"/*, "Phobos" */) );

    dependentVariablesList.push_back(
                std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "AlienSpaceship", "Phobos" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, alienSpaceshipInitialState, simulationEndEpoch,
                cowell, dependentVariablesToSave );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "AlienSpaceship", alienSpaceshipInitialState, "Phobos" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "AlienSpaceship", radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
                dependentVariablesToSave );

    // Retrieve dependent variables history.
    std::map< double, Eigen::VectorXd > dependentVariablesHistory = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getDependentVariableHistory( );

    // Retrieve dependent variables interface.
    std::shared_ptr< SingleArcDependentVariablesInterface > dependentVariablesInterface =
            std::dynamic_pointer_cast< SingleArcDependentVariablesInterface >(
            variationalEquationsSimulator.getDependentVariablesInterface( ) );

    // Create dependent variables interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory ), 4 );

    std::vector< double > testEpochs;
    testEpochs.push_back( simulationStartEpoch );
    testEpochs.push_back( ( simulationEndEpoch -  simulationStartEpoch ) / 4.0 );
    testEpochs.push_back( ( simulationEndEpoch -  simulationStartEpoch ) / 2.0 );
    testEpochs.push_back( 3.0 * ( simulationEndEpoch -  simulationStartEpoch ) / 4.0 );
    testEpochs.push_back( simulationEndEpoch );

    // Check consistency between interpolator results and interface results.
    for ( unsigned int i = 0 ; i < testEpochs.size( ) ; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariablesInterpolator->interpolate( testEpochs[ i ] ),
                                           dependentVariablesInterface->getDependentVariables( testEpochs[ i ] ),
                                           std::numeric_limits< double >::epsilon( ) );
    }


    std::map< double, Eigen::VectorXd > totalAccelerationHistory;
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ; itr != dependentVariablesHistory.end( ) ; itr++ )
    {
        totalAccelerationHistory[ itr->first ] = itr->second.segment( 0, 3 );
    }

    // Create total acceleration history interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

    // Total acceleration dependent variable settings.
    std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
            = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );

    // Check consistency between interpolator results and interface results, for a single dependent variable.
    for ( unsigned int i = 0 ; i < testEpochs.size( ) ; i++ )
    {
//        std::cout << "TEST " << i << "\n\n";
//        std::cout << "interpolated dependent variable: " << totalAccelerationInterpolator->interpolate( testEpochs[ i ] ).transpose( ) << "\n\n";
//        std::cout << "dependent variable from interface: " << dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ i ] ).transpose( ) << "\n\n";
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ i ] ),
                                           dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ i ] ),
                                           std::numeric_limits< double >::epsilon( ) );
    }

}



//! Test the dependent variables computation against their numerical propagation.
/*!
 *  Test the dependent variables computation against their numerical propagation. This unit test
 *  ...
 */
BOOST_AUTO_TEST_CASE( testMultiArcDependentVariablesInterface )
{

//    // Load Spice kernels.
//    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 8.0 * tudat::physical_constants::JULIAN_DAY;


    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );


    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "ECLIPJ2000" );
    }
    bodySettings[ "Phobos" ] = std::make_shared< BodySettings >( );
    bodySettings[ "Phobos" ]->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
    bodySettings[ "Phobos" ]->gravityFieldSettings =
            getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    NamedBodyMap bodyMap = createBodies( bodySettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "AlienSpaceship" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "AlienSpaceship" ]->setConstantBodyMass( 400.0 );

//    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
//                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                            < double, Eigen::Vector6d > >( ), "Phobos", "ECLIPJ2000" ) );

    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ), "Phobos", "ECLIPJ2000" ) );

//    // Create aerodynamic coefficient interface settings.
//    double referenceArea = 4.0;
//    double aerodynamicCoefficient = 1.2;
//    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
//            std::make_shared< ConstantAerodynamicCoefficientSettings >(
//                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

//    // Create and set aerodynamic coefficients object
//    bodyMap[ "AlienSpaceship" ]->setAerodynamicCoefficientInterface(
//                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "AlienSpaceship" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > alienSpaceshipRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "AlienSpaceship" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    alienSpaceshipRadiationPressureSettings, "AlienSpaceship", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAlienSpaceship;
    accelerationsOfAlienSpaceship[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAlienSpaceship[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAlienSpaceship[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
//    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationMap[ "AlienSpaceship" ] = accelerationsOfAlienSpaceship;
    bodiesToPropagate.push_back( "AlienSpaceship" );
    centralBodies.push_back( "Phobos" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d alienSpaceshipInitialStateInKeplerianElements;
    alienSpaceshipInitialStateInKeplerianElements( semiMajorAxisIndex ) = 100.0E3;
    alienSpaceshipInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
    alienSpaceshipInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
    alienSpaceshipInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    alienSpaceshipInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    alienSpaceshipInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    double phobosGravitationalParameter = bodyMap.at( "Phobos" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Define arc length
    double arcDuration = 3.0 * physical_constants::JULIAN_DAY;
//    double arcOverlap = 3600.0;


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;


    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_acceleration_dependent_variable, "AlienSpaceship" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Sun", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Mars", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Phobos", point_mass_gravity, "AlienSpaceship", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "Earth", third_body_spherical_harmonic_gravity, "AlienSpaceship", "Phobos" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
//                    "AlienSpaceship", "Sun", cannon_ball_radiation_pressure, "AlienSpaceship"/*, "Phobos" */) );

    dependentVariablesList.push_back(
                std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >(
                    "AlienSpaceship", "AlienSpaceship", "Phobos" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


    // Create propagator settings (including initial state taken from Kepler orbit) for each arc
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    std::vector< double > arcStartTimes;
    std::vector< double > arcEndTimes;
    double currentTime = simulationStartEpoch;
    while( currentTime <= simulationEndEpoch )
    {
        arcStartTimes.push_back( currentTime );

        Eigen::Vector6d currentArcInitialState = convertKeplerianToCartesianElements(
                    propagateKeplerOrbit( alienSpaceshipInitialStateInKeplerianElements, currentTime - simulationStartEpoch,
                                          phobosGravitationalParameter ), phobosGravitationalParameter );
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                                              centralBodies, accelerationModelMap, bodiesToPropagate, currentArcInitialState,
                                              currentTime + arcDuration, cowell, dependentVariablesToSave ) );

        arcEndTimes.push_back( currentTime + arcDuration );

        currentTime += arcDuration;
    }

    // Create propagator settings
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );

//    const Eigen::Vector6d alienSpaceshipInitialState = convertKeplerianToCartesianElements(
//                alienSpaceshipInitialStateInKeplerianElements, phobosGravitationalParameter );





//    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
//            std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                centralBodies, accelerationModelMap, bodiesToPropagate, alienSpaceshipInitialState, simulationEndEpoch,
//                cowell, dependentVariablesToSave );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create concatenated list of arc initial states
    Eigen::VectorXd alienSpaceshipInitialState = Eigen::VectorXd( 6 * arcStartTimes.size( ) );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        alienSpaceshipInitialState.segment( i * 6, 6 ) = propagatorSettingsList.at( i )->getInitialStates( );
    }

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "AlienSpaceship", alienSpaceshipInitialState, arcStartTimes, "Phobos" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "AlienSpaceship", radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create simulation object and propagate dynamics.
//    SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
//                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
//                dependentVariablesToSave );

    // Create dynamics simulator
    MultiArcVariationalEquationsSolver< > variationalEquationsSimulator =
            MultiArcVariationalEquationsSolver< >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, arcStartTimes, true, std::shared_ptr< IntegratorSettings< double > >( ), false, true, true,
                dependentVariablesToSave );

    // Retrieve dependent variables history.
    std::vector< std::map< double, Eigen::VectorXd > > dependentVariablesHistory = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getDependentVariableHistory( );

    // Retrieve dependent variables interface.
    std::shared_ptr< MultiArcDependentVariablesInterface > dependentVariablesInterface =
            std::dynamic_pointer_cast< MultiArcDependentVariablesInterface >(
            variationalEquationsSimulator.getDependentVariablesInterface( ) );

    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        // Create dependent variables interpolator.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory[ i ] ),
                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory[ i ] ), 4 );

        std::vector< double > testEpochs;
        testEpochs.push_back( arcStartTimes[ i ] + 10.0 );
        testEpochs.push_back( arcStartTimes[ i ] + ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 4.0 );
        testEpochs.push_back( arcStartTimes[ i ] + ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 2.0 );
        testEpochs.push_back( arcStartTimes[ i ] + 3.0 * ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 4.0 );
        testEpochs.push_back( arcEndTimes[ i ] - 10.0 );

        // Check consistency between interpolator results and interface results.
        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
        {
//            std::cout << "interpolator: " << dependentVariablesInterpolator->interpolate( testEpochs[ j ] ).transpose( ) << "\n\n";
//            std::cout << "interface: " << dependentVariablesInterface->getDependentVariables( testEpochs[ j ] ).transpose( ) << "\n\n";
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariablesInterpolator->interpolate( testEpochs[ i ] ),
                                               dependentVariablesInterface->getDependentVariables( testEpochs[ i ] ),
                                               std::numeric_limits< double >::epsilon( ) );
        }


        std::map< double, Eigen::VectorXd > totalAccelerationHistory;
        for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory[ i ].begin( ) ; itr != dependentVariablesHistory[ i ].end( ) ; itr++ )
        {
            totalAccelerationHistory[ itr->first ] = itr->second.segment( 0, 3 );
        }


        // Create total acceleration history interpolator.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

        // Total acceleration dependent variable settings.
        std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
                = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );

        // Partial of total acceleration w.r.t. translational state dependent variable settings.
        std::shared_ptr< TotalAccelerationPartialWrtStateSaveSettings > partialTotalAccelerationWrtStateDependentVariable
                = std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "AlienSpaceship", "AlienSpaceship", "Phobos" );

        // Check consistency between interpolator results and interface results, for a single dependent variable.
        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
        {
//            std::cout << "TEST " << j << "\n\n";
//            std::cout << "interpolated dependent variable: " << totalAccelerationInterpolator->interpolate( testEpochs[ j ] ).transpose( ) << "\n\n";
//            std::cout << "dependent variable from interface: " << dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ).transpose( ) << "\n\n";
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ j ] ),
                                               dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ),
                                               std::numeric_limits< double >::epsilon( ) );

            std::cout << "partial of total acceleration w.r.t. state from interface: " <<
                         dependentVariablesInterface->getSingleDependentVariable( partialTotalAccelerationWrtStateDependentVariable, testEpochs[ j ] ).transpose( ) << "\n\n";
        }


    }

}


////! Test the dependent variables computation against their numerical propagation.
///*!
// *  Test the dependent variables computation against their numerical propagation. This unit test
// *  ...
// */
//BOOST_AUTO_TEST_CASE( testHybridArcDependentVariablesInterface )
//{

////    // Load Spice kernels.
////    spice_interface::loadStandardSpiceKernels( );

//    // Set simulation time settings.
//    const double simulationStartEpoch = 0.0;
//    const double simulationEndEpoch = 8.0 * tudat::physical_constants::JULIAN_DAY;


//    // Load spice kernel.
//    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );


//    // Define body settings for simulation.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );
//    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Mars" );

//    // Create body objects.
//    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
//            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
//    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
//    {
//        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
//        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "ECLIPJ2000" );
//    }
//    bodySettings[ "Phobos" ] = std::make_shared< BodySettings >( );
//    bodySettings[ "Phobos" ]->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
//    bodySettings[ "Phobos" ]->gravityFieldSettings =
//            getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
//    NamedBodyMap bodyMap = createBodies( bodySettings );



//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////             DEFINE TABULATED PHOBOS EPHEMERIS               ///////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////    Eigen::Vector6d initialStatePhobos = bodyMap[ "Phobos" ]->getEphemeris( )->getCartesianState( simulationStartEpoch );

////    bodyMap[ "Phobos" ]->setEphemeris(  std::make_shared< TabulatedCartesianEphemeris< > >(
////                                            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Mars", "ECLIPJ2000" ) );

////    std::shared_ptr< Ephemeris > phobosEphemeris =
////            createBodyEphemeris( getDefaultEphemerisSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 ), "Phobos" );

////    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > defaultPhobosStateInterpolator =
////            std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( phobosEphemeris )->getInterpolator( );

////    std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
////                bodyMap.at( "Phobos" )->getEphemeris( ) )->resetInterpolator( defaultPhobosStateInterpolator );



//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create spacecraft object.
//    bodyMap[ "AlienSpaceship" ] = std::make_shared< simulation_setup::Body >( );
//    bodyMap[ "AlienSpaceship" ]->setConstantBodyMass( 400.0 );

////    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
////                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
////                                            < double, Eigen::Vector6d > >( ), "Phobos", "ECLIPJ2000" ) );

//    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
//                                            std::map< double, std::shared_ptr< Ephemeris > >( ), "Phobos", "ECLIPJ2000" ) );

////    // Create aerodynamic coefficient interface settings.
////    double referenceArea = 4.0;
////    double aerodynamicCoefficient = 1.2;
////    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
////            std::make_shared< ConstantAerodynamicCoefficientSettings >(
////                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

////    // Create and set aerodynamic coefficients object
////    bodyMap[ "AlienSpaceship" ]->setAerodynamicCoefficientInterface(
////                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "AlienSpaceship" ) );

//    // Create radiation pressure settings
//    double referenceAreaRadiation = 4.0;
//    double radiationPressureCoefficient = 1.2;
//    std::vector< std::string > occultingBodies;
//    occultingBodies.push_back( "Mars" );
//    std::shared_ptr< RadiationPressureInterfaceSettings > alienSpaceshipRadiationPressureSettings =
//            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
//                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

//    // Create and set radiation pressure settings
//    bodyMap[ "AlienSpaceship" ]->setRadiationPressureInterface(
//                "Sun", createRadiationPressureInterface(
//                    alienSpaceshipRadiationPressureSettings, "AlienSpaceship", bodyMap ) );


//    // Finalize body creation.
//    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );



////    // Tabulated states derived from default ephemeris
////    std::map< double, Eigen::Vector6d > tabulatedStatesPhobos;
////    double currentTimeTabulatedStates = simulationStartEpoch;
////    while( currentTimeTabulatedStates <= simulationEndEpoch )
////    {
////        tabulatedStatesPhobos[ currentTimeTabulatedStates ] = bodyMap[ "Phobos" ]->getEphemeris( )->getCartesianState( currentTimeTabulatedStates ); // convertKeplerianToCartesianElements(
//////                    propagateKeplerOrbit( keplerianElementsPrimary, currentTimeTabulatedStates - initialEphemerisTime,
//////                                          gravitationalParameterSun ), gravitationalParameterSun );
////        currentTimeTabulatedStates += 10.0;
////    }

////    bodyMap[ "Phobos" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
////                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Vector6d > >( tabulatedStatesPhobos, 6 ), "SSB", "ECLIPJ2000" ) );
//////    bodySettings[ "Phobos" ]->ephemerisSettings = std::make_shared< TabulatedEphemerisSettings >(
//////                tabulatedStatesPhobos, "SSB", "ECLIPJ2000" );



//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Set accelerations for Phobos
//    SelectedAccelerationMap singleArcAccelerationMap;
//    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
//    accelerationsOfPhobos[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    accelerationsOfPhobos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
//    singleArcAccelerationMap[ "Phobos" ] = accelerationsOfPhobos;

//    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
//    singleArcBodiesToIntegrate.push_back( "Phobos" );
//    singleArcCentralBodies.push_back( "SSB" );

//    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
//                bodyMap, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );


//    // Define propagator settings variables.
//    SelectedAccelerationMap accelerationMap;
//    std::vector< std::string > bodiesToPropagate;
//    std::vector< std::string > centralBodies;

//    // Define propagation settings.
//    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAlienSpaceship;
//    accelerationsOfAlienSpaceship[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

//    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                   basic_astrodynamics::central_gravity ) );
//    accelerationsOfAlienSpaceship[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     basic_astrodynamics::central_gravity ) );
//    accelerationsOfAlienSpaceship[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     basic_astrodynamics::central_gravity ) );
//    accelerationsOfAlienSpaceship[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

//    accelerationMap[ "AlienSpaceship" ] = accelerationsOfAlienSpaceship;
//    bodiesToPropagate.push_back( "AlienSpaceship" );
//    centralBodies.push_back( "Phobos" );

//    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Set Keplerian elements for Asterix.
//    Eigen::Vector6d alienSpaceshipInitialStateInKeplerianElements;
//    alienSpaceshipInitialStateInKeplerianElements( semiMajorAxisIndex ) = 100.0E3;
//    alienSpaceshipInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
//    alienSpaceshipInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
//    alienSpaceshipInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
//    alienSpaceshipInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
//    alienSpaceshipInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

//    double phobosGravitationalParameter = bodyMap.at( "Phobos" )->getGravityFieldModel( )->getGravitationalParameter( );

//    // Define arc length
//    double arcDuration = 3.0 * physical_constants::JULIAN_DAY;
////    double arcOverlap = 3600.0;


//    // Define list of dependent variables.
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >
//                                      ( total_acceleration_dependent_variable, "Phobos" ) );
//    dependentVariablesList.push_back( std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "Phobos", "Phobos" ) );
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
//                    total_acceleration_dependent_variable, "AlienSpaceship" ) );
//    dependentVariablesList.push_back( std::make_shared< AccelerationPartialWrtStateSaveSettings >(
//                    "AlienSpaceship", "Sun", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );

//    // Define list of dependent variables to save for single-arc dynamics simulator.
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > singleArcDependentVariablesList;
//    singleArcDependentVariablesList.push_back( dependentVariablesList[ 0 ] );
//    singleArcDependentVariablesList.push_back( dependentVariablesList[ 1 ] );

////    singleArcDependentVariablesList.push_back(
////                std::make_shared< SingleDependentVariableSaveSettings >(
////                    total_acceleration_dependent_variable, "Phobos" ) );
////    singleArcDependentVariablesList.push_back(
////                std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >(
////                    "Phobos", "Phobos" ) );


//    // Define list of dependent variables to save for multi-arc dynamics simulator.
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > multiArcDependentVariablesList;
//    multiArcDependentVariablesList.push_back( dependentVariablesList[ 2 ] );
//    multiArcDependentVariablesList.push_back( dependentVariablesList[ 3 ] );

////    dependentVariablesList.push_back(
////                std::make_shared< SingleDependentVariableSaveSettings >(
////                    total_acceleration_dependent_variable, "AlienSpaceship" ) );
////    dependentVariablesList.push_back(
////                std::make_shared< AccelerationPartialWrtStateSaveSettings >(
////                    "AlienSpaceship", "Sun", third_body_central_gravity, "AlienSpaceship", "Phobos" ) );

//    // Create object with list of dependent variables
//    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
//            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


//    // Create single-arc propagator settings for Phobos
//    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
//                singleArcBodiesToIntegrate, singleArcCentralBodies, bodyMap, simulationStartEpoch );
//    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
//            std::make_shared< TranslationalStatePropagatorSettings< > >(
//                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
//                singleArcInitialStates, simulationEndEpoch, cowell, std::make_shared< DependentVariableSaveSettings >( singleArcDependentVariablesList ) );


//    // Create propagator settings (including initial state taken from Kepler orbit) for each arc
//    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//    std::vector< double > arcStartTimes;
//    std::vector< double > arcEndTimes;
//    double currentTime = simulationStartEpoch;
//    while( currentTime <= simulationEndEpoch )
//    {
//        arcStartTimes.push_back( currentTime );

//        Eigen::Vector6d currentArcInitialState = convertKeplerianToCartesianElements(
//                    propagateKeplerOrbit( alienSpaceshipInitialStateInKeplerianElements, currentTime - simulationStartEpoch,
//                                          phobosGravitationalParameter ), phobosGravitationalParameter );
//        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                                              centralBodies, accelerationModelMap, bodiesToPropagate, currentArcInitialState,
//                                              currentTime + arcDuration, cowell, std::make_shared< DependentVariableSaveSettings >( multiArcDependentVariablesList ) ) );

//        arcEndTimes.push_back( currentTime + arcDuration );

//        currentTime += arcDuration;
//    }

//    // Create multi-arc propagator settings
//    std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
//            std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );


//    // Create hybrid propagator settings.
//    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
//                singleArcPropagatorSettings, multiArcPropagatorSettings );



////    const Eigen::Vector6d alienSpaceshipInitialState = convertKeplerianToCartesianElements(
////                alienSpaceshipInitialStateInKeplerianElements, phobosGravitationalParameter );





////    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
////            std::make_shared< TranslationalStatePropagatorSettings< double > >(
////                centralBodies, accelerationModelMap, bodiesToPropagate, alienSpaceshipInitialState, simulationEndEpoch,
////                cowell, dependentVariablesToSave );

//    const double fixedStepSize = 10.0;
//    std::shared_ptr< IntegratorSettings< > > integratorSettings =
//            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, fixedStepSize );

//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create concatenated list of arc initial states
//    Eigen::VectorXd alienSpaceshipInitialState = Eigen::VectorXd( 6 * arcStartTimes.size( ) );
//    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
//    {
//        alienSpaceshipInitialState.segment( i * 6, 6 ) = propagatorSettingsList.at( i )->getInitialStates( );
//    }

//    // Define list of parameters to estimate.
//    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//    parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
//                                  "AlienSpaceship", alienSpaceshipInitialState, arcStartTimes, "Phobos" ) );
//    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                    singleArcBodiesToIntegrate.at( 0 ), singleArcInitialStates, singleArcCentralBodies.at( 0 ) ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "AlienSpaceship", radiation_pressure_coefficient ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
//    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
//    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

//    // Create parameters
//    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//            createParametersToEstimate( parameterNames, bodyMap );

//    // Print identifiers and indices of parameters to terminal.
//    printEstimatableParameterEntries( parametersToEstimate );




////    // Load spice kernels.
////    spice_interface::loadStandardSpiceKernels( );

////    // Create list of bodies to create.
////    std::vector< std::string > bodyNames;
////    bodyNames.push_back( "Jupiter" );
////    bodyNames.push_back( "Mars" );
//////    bodyNames.push_back( "Europa" );
//////    bodyNames.push_back( "Ganymede" );

////    // Specify initial time
////    double initialTime = 0.0;
////    double finalTime = 4.0 * 86400.0;

////    // Get body settings.
////    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
////            getDefaultBodySettings( bodyNames, initialTime - 3600.0, finalTime + 3600.0 );

////    // Create bodies needed in simulation
////    NamedBodyMap bodyMap = createBodies( bodySettings );
////    bodyMap[ "Spacecraft" ] = std::make_shared< Body >( );
////    bodyMap[ "Spacecraft" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
////                                               std::map< double, std::shared_ptr< Ephemeris > >( ), "Jupiter", "ECLIPJ2000" ) );
////    setGlobalFrameBodyEphemerides( bodyMap, "Jupiter", "ECLIPJ2000" );

////    SelectedAccelerationMap singleArcAccelerationMap;
////    std::vector< std::string > singleArcBodiesToPropagate = { "Mars" }; //, "Europa", "Ganymede" };
////    std::vector< std::string > singleArcCentralBodies = { "Sun" }; //, "Jupiter", "Jupiter" };

////    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
////    {
////        singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ "Sun" ].push_back(
////                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

//////        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
//////        {
//////            if( i != j )
//////            {
//////                singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ singleArcBodiesToPropagate.at( j ) ].push_back(
//////                            std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
//////            }
//////        }
////    }

////    basic_astrodynamics::AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
////                bodyMap, singleArcAccelerationMap, singleArcBodiesToPropagate, singleArcCentralBodies );

////    Eigen::VectorXd singleArcInitialState = getInitialStatesOfBodies(
////                singleArcBodiesToPropagate, singleArcCentralBodies, bodyMap, initialTime );

////    std::shared_ptr< TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings =
////            std::make_shared< TranslationalStatePropagatorSettings< double > >
////            ( singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToPropagate,singleArcInitialState, finalTime );

////    std::vector< std::string > multiArcBodiesToPropagate =
////    { "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft" };

////    std::vector< std::string > multiArcCentralBodies =
////    { "Io", "Ganymede", "Europa", "Ganymede", "Io", "Europa" };

////    std::vector< double > arcStartTimes =
////    { 3600.0, 3.0 * 3600, 5.0 * 3600.0, 7.0 * 3600.0, 9.0 * 3600.0, 11.0 * 3600.0 };
////    std::map< std::string, std::vector< double > > arcStartTimesPerBody;

////    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
////    {
////        arcStartTimesPerBody[ multiArcCentralBodies.at( i ) ].push_back( arcStartTimes.at( i ) );
////    }
////    double arcDuration = 3600.0;


////    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
////    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > multiArcPropagationSettingsList;
////    std::map< std::string, std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > >
////            multiArcPropagationSettingsListPerCentralBody;
////    std::map< std::string, std::vector< int > > perBodyIndicesInFullPropagation;

////    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
////    {
////        Eigen::Vector6d spacecraftInitialStateInKeplerianElements;
////        spacecraftInitialStateInKeplerianElements( semiMajorAxisIndex ) = 3500.0E3;
////        spacecraftInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
////        spacecraftInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians(
////                    static_cast< double >( i * 30 ) );
////        spacecraftInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians(
////                    static_cast< double >( i * 30 ) );
////        spacecraftInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians(
////                    static_cast< double >( i * 30 ) );
////        spacecraftInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians(
////                    static_cast< double >( i * 30 ) );
////        double centralBodyGravitationalParameter = bodyMap.at( multiArcCentralBodies.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );
////        multiArcSystemInitialStates.push_back( convertKeplerianToCartesianElements(
////                                                   spacecraftInitialStateInKeplerianElements, centralBodyGravitationalParameter ) );

////        SelectedAccelerationMap multiArcAccelerationMap;

////        multiArcAccelerationMap[ "Spacecraft" ][ "Jupiter" ].push_back(
////                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

////        for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
////        {
////            multiArcAccelerationMap[ "Spacecraft" ][ singleArcBodiesToPropagate.at( i ) ].push_back(
////                        std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
////        }

////        basic_astrodynamics::AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
////                    bodyMap, multiArcAccelerationMap, { multiArcBodiesToPropagate.at( i ) }, { multiArcCentralBodies.at( i ) } );

////        multiArcPropagationSettingsList.push_back(
////                    std::make_shared< TranslationalStatePropagatorSettings< double > >
////                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
////                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
////                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration ) );
////        multiArcPropagationSettingsListPerCentralBody[
////                multiArcCentralBodies.at( i ) ].push_back(
////                    std::make_shared< TranslationalStatePropagatorSettings< double > >
////                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
////                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
////                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration ) );
////        perBodyIndicesInFullPropagation[ multiArcCentralBodies.at( i ) ].push_back( i );

////    }
////    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagationSettings =
////            std::make_shared< MultiArcPropagatorSettings< > >( multiArcPropagationSettingsList );

////    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
////            std::make_shared< HybridArcPropagatorSettings< > >( singleArcPropagatorSettings, multiArcPropagationSettings );

////    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
////    parameterNames.push_back(
////                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
////                    multiArcBodiesToPropagate.at( 0 ), multiArcPropagationSettings->getInitialStates( ),
////                    arcStartTimes, multiArcCentralBodies ) );
////    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
////    {
////        parameterNames.push_back(
////                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
////                        singleArcBodiesToPropagate.at( i ), singleArcInitialState.segment( 6 * i, 6 ), singleArcCentralBodies.at( i ) ) );
////        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
////                                      singleArcBodiesToPropagate.at( i ), gravitational_parameter ) );
////    }

////    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
////            createParametersToEstimate< double >( parameterNames, bodyMap, hybridArcPropagatorSettings );
////    printEstimatableParameterEntries( parametersToEstimate );

////    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings =
////            std::make_shared< IntegratorSettings< > >
////            ( rungeKutta4, initialTime, 60.0 );

////    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings =
////            std::make_shared< IntegratorSettings< > >
////            ( rungeKutta4, TUDAT_NAN, 15.0 );

////    // Create dynamics simulator
////    HybridArcVariationalEquationsSolver< > variationalEquations =
////            HybridArcVariationalEquationsSolver< >(
////                bodyMap, singleArcIntegratorSettings, multiArcIntegratorSettings,
////                hybridArcPropagatorSettings, parametersToEstimate, arcStartTimes, true, false, true );



//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////    // Create simulation object and propagate dynamics.
////    SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
////                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
////                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
////                dependentVariablesToSave );

////    // Create dynamics simulator
////    MultiArcVariationalEquationsSolver< > variationalEquationsSimulator =
////            MultiArcVariationalEquationsSolver< >(
////                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, arcStartTimes, true, std::shared_ptr< IntegratorSettings< double > >( ), false, true, true,
////                dependentVariablesToSave );

////    // Create hybrid arc variational equations solver.
////    HybridArcVariationalEquationsSolver< > variationalEquationsSimulator( bodyMap, integratorSettings, hybridArcPropagatorSettings, parametersToEstimate, arcStartTimes,
////                                                                          true, false, true, dependentVariablesToSave );


//    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings =
//            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, 10.0 );

//    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings =
//            std::make_shared< IntegratorSettings< > >( rungeKutta4, TUDAT_NAN, 10.0 );

////    // Create dynamics simulator
////    HybridArcVariationalEquationsSolver< > variationalEquations =
////            HybridArcVariationalEquationsSolver< >(
////                bodyMap, singleArcIntegratorSettings, multiArcIntegratorSettings,
////                hybridArcPropagatorSettings, parametersToEstimate, arcStartTimes, true, false, true );

////    // Retrieve dependent variables history.
////    std::vector< std::map< double, Eigen::VectorXd > > multiArcDependentVariablesHistory = variationalEquationsSimulator.getMultiArcSolver( )->getDynamicsSimulator( )
////            ->getDependentVariableHistory( );
////    std::map< double, Eigen::VectorXd > singleArcDependentVariablesHistory = variationalEquationsSimulator.getSingleArcSolver( )->getDynamicsSimulator( )->getDependentVariableHistory( );

////    // Retrieve dependent variables interface.
////    std::shared_ptr< HybridArcDependentVariablesInterface > dependentVariablesInterface = std::dynamic_pointer_cast< HybridArcDependentVariablesInterface >(
////            variationalEquationsSimulator.getDependentVariablesInterface( ) );

////    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
////    {
////        // Create dependent variables interpolator.
////        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
////                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
////                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( multiArcDependentVariablesHistory[ i ] ),
////                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( multiArcDependentVariablesHistory[ i ] ), 4 );

////        std::vector< double > testEpochs;
////        testEpochs.push_back( arcStartTimes[ i ] + 10.0 );
////        testEpochs.push_back( arcStartTimes[ i ] + ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 4.0 );
////        testEpochs.push_back( arcStartTimes[ i ] + ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 2.0 );
////        testEpochs.push_back( arcStartTimes[ i ] + 3.0 * ( arcEndTimes[ i ] -  arcStartTimes[ i ] ) / 4.0 );
////        testEpochs.push_back( arcEndTimes[ i ] - 10.0 );

////        // Check consistency between interpolator results and interface results.
////        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
////        {
////            std::cout << "interpolator: " << dependentVariablesInterpolator->interpolate( testEpochs[ j ] ).transpose( ) << "\n\n";
////            std::cout << "interface: " << dependentVariablesInterface->getDependentVariables( testEpochs[ j ] ).transpose( ) << "\n\n";
//////            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariablesInterpolator->interpolate( testEpochs[ i ] ),
//////                                               dependentVariablesInterface->getDependentVariables( testEpochs[ i ] ),
//////                                               std::numeric_limits< double >::epsilon( ) );
////        }


////        std::map< double, Eigen::VectorXd > totalAccelerationHistory;
////        for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory[ i ].begin( ) ; itr != dependentVariablesHistory[ i ].end( ) ; itr++ )
////        {
////            totalAccelerationHistory[ itr->first ] = itr->second.segment( 0, 3 );
////        }


////        // Create total acceleration history interpolator.
////        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
////                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
////                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
////                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

////        // Total acceleration dependent variable settings.
////        std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
////                = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );

////        // Check consistency between interpolator results and interface results, for a single dependent variable.
////        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
////        {
//////            std::cout << "TEST " << j << "\n\n";
//////            std::cout << "interpolated dependent variable: " << totalAccelerationInterpolator->interpolate( testEpochs[ j ] ).transpose( ) << "\n\n";
//////            std::cout << "dependent variable from interface: " << dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ).transpose( ) << "\n\n";
////            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ j ] ),
////                                               dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ),
////                                               std::numeric_limits< double >::epsilon( ) );
////        }


////    }

//}


BOOST_AUTO_TEST_CASE( testHybridArcDependentVariablesInterface2 )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Phobos" );
//    bodyNames.push_back( "Europa" );
//    bodyNames.push_back( "Ganymede" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 4.0 * 86400.0;

    // Get body settings.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialTime - 300.0, finalTime + 300.0 );

    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ "AlienSpaceship" ] = std::make_shared< Body >( );
    bodyMap[ "AlienSpaceship" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                               std::map< double, std::shared_ptr< Ephemeris > >( ), "Mars", "ECLIPJ2000" ) );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    SelectedAccelerationMap singleArcAccelerationMap;
    std::vector< std::string > singleArcBodiesToPropagate = { "Phobos" }; //, "Europa", "Ganymede" };
    std::vector< std::string > singleArcCentralBodies = { "Mars" }; //, "Jupiter", "Jupiter" };

    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
    {
        singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ "Mars" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
        {
            if( i != j )
            {
                singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ singleArcBodiesToPropagate.at( j ) ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
            }
        }
    }

    basic_astrodynamics::AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, singleArcAccelerationMap, singleArcBodiesToPropagate, singleArcCentralBodies );

    Eigen::VectorXd singleArcInitialState = getInitialStatesOfBodies(
                singleArcBodiesToPropagate, singleArcCentralBodies, bodyMap, initialTime );

//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >
//                                      ( total_acceleration_dependent_variable, "Phobos" ) );
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >
//                                      ( total_acceleration_dependent_variable, "Phobos" ) );

    // Define list of dependent variables.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >
                                      ( total_acceleration_dependent_variable, "Phobos" ) );
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >
//                                      ( relative_position_dependent_variable, "Phobos", "Mars" ) );
//    dependentVariablesList.push_back( std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "Phobos", "Phobos", "Mars" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                    total_acceleration_dependent_variable, "AlienSpaceship" ) );
//    dependentVariablesList.push_back( std::make_shared< AccelerationPartialWrtStateSaveSettings >(
//                    "AlienSpaceship", "Phobos", point_mass_gravity, "AlienSpaceship"/*, "Phobos"*/ ) );

    // Define list of dependent variables to save for single-arc dynamics simulator.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > singleArcDependentVariablesList;
    singleArcDependentVariablesList.push_back( dependentVariablesList[ 0 ] );
//    singleArcDependentVariablesList.push_back( dependentVariablesList[ 1 ] );

    // Define list of dependent variables to save for multi-arc dynamics simulator.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > multiArcDependentVariablesList;
    multiArcDependentVariablesList.push_back( dependentVariablesList[ 1 ] );
//    multiArcDependentVariablesList.push_back( dependentVariablesList[ 2 ] );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

//    // Create object with list of dependent variables
//    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
//            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToPropagate,singleArcInitialState, finalTime,
              cowell, std::make_shared< DependentVariableSaveSettings >( singleArcDependentVariablesList ) );

    std::vector< std::string > multiArcBodiesToPropagate = { "AlienSpaceship", "AlienSpaceship", "AlienSpaceship", "AlienSpaceship", "AlienSpaceship", "AlienSpaceship" };
    std::vector< std::string > multiArcCentralBodies = { "Phobos", "Phobos", "Phobos", "Phobos", "Phobos", "Phobos" };

    std::vector< double > arcStartTimes =
    { 3600.0, 3.0 * 3600, 5.0 * 3600.0, 7.0 * 3600.0, 9.0 * 3600.0, 11.0 * 3600.0 };
    std::map< std::string, std::vector< double > > arcStartTimesPerBody;

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        arcStartTimesPerBody[ multiArcCentralBodies.at( i ) ].push_back( arcStartTimes.at( i ) );
    }
    double arcDuration = 3600.0;


    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > multiArcPropagationSettingsList;
    std::map< std::string, std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > >
            multiArcPropagationSettingsListPerCentralBody;
    std::map< std::string, std::vector< int > > perBodyIndicesInFullPropagation;

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        Eigen::Vector6d spacecraftInitialStateInKeplerianElements;
        spacecraftInitialStateInKeplerianElements( semiMajorAxisIndex ) = 3500.0E3;
        spacecraftInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        spacecraftInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        double centralBodyGravitationalParameter = bodyMap.at( multiArcCentralBodies.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );
        multiArcSystemInitialStates.push_back( convertKeplerianToCartesianElements(
                                                   spacecraftInitialStateInKeplerianElements, centralBodyGravitationalParameter ) );

        SelectedAccelerationMap multiArcAccelerationMap;

        multiArcAccelerationMap[ "AlienSpaceship" ][ "Mars" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

        for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
        {
            multiArcAccelerationMap[ "AlienSpaceship" ][ singleArcBodiesToPropagate.at( i ) ].push_back(
                        std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
        }

        basic_astrodynamics::AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                    bodyMap, multiArcAccelerationMap, { multiArcBodiesToPropagate.at( i ) }, { multiArcCentralBodies.at( i ) } );

        multiArcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration, cowell,
                      std::make_shared< DependentVariableSaveSettings >( multiArcDependentVariablesList ) ) );
        multiArcPropagationSettingsListPerCentralBody[
                multiArcCentralBodies.at( i ) ].push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration, cowell,
                      std::make_shared< DependentVariableSaveSettings >( multiArcDependentVariablesList ) ) );
        perBodyIndicesInFullPropagation[ multiArcCentralBodies.at( i ) ].push_back( i );

    }
    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagationSettings =
            std::make_shared< MultiArcPropagatorSettings< > >( multiArcPropagationSettingsList );

    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
            std::make_shared< HybridArcPropagatorSettings< > >( singleArcPropagatorSettings, multiArcPropagationSettings );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    multiArcBodiesToPropagate.at( 0 ), multiArcPropagationSettings->getInitialStates( ),
                    arcStartTimes, multiArcCentralBodies ) );
    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        singleArcBodiesToPropagate.at( i ), singleArcInitialState.segment( 6 * i, 6 ), singleArcCentralBodies.at( i ) ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                                      singleArcBodiesToPropagate.at( i ), gravitational_parameter ) );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodyMap, hybridArcPropagatorSettings );
    printEstimatableParameterEntries( parametersToEstimate );

    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialTime, 60.0 );

    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, TUDAT_NAN, 15.0 );

    // Create dynamics simulator
    HybridArcVariationalEquationsSolver< > variationalEquations =
            HybridArcVariationalEquationsSolver< >(
                bodyMap, singleArcIntegratorSettings, multiArcIntegratorSettings,
                hybridArcPropagatorSettings, parametersToEstimate, arcStartTimes, true, false, true, dependentVariablesToSave );

    // Retrieve dependent variables interface.
    std::shared_ptr< HybridArcDependentVariablesInterface > dependentVariablesInterface =
            std::dynamic_pointer_cast< HybridArcDependentVariablesInterface >(
            variationalEquations.getDependentVariablesInterface( ) );

    std::cout << "dependent variables from interface: " << dependentVariablesInterface->getDependentVariables( ( finalTime - initialTime ) / 2.0 ).transpose( ) << "\n\n";

//    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > fullMultiArcVariationalSolution =
//            variationalEquations.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
//    std::vector< std::map< double, Eigen::VectorXd > > fullMultiArcStateSolution =
//            variationalEquations.getMultiArcSolver( )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

//    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
//    {
//        std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPerBodyPropagationSettings =
//                std::make_shared< MultiArcPropagatorSettings< > >( multiArcPropagationSettingsListPerCentralBody.at(
//                                                                       singleArcBodiesToPropagate.at( i ) ) );
//        std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPerBodyPropagatorSettings =
//                std::make_shared< HybridArcPropagatorSettings< > >(
//                    singleArcPropagatorSettings, multiArcPerBodyPropagationSettings );

//        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNamesPerBody;
//        parameterNamesPerBody.push_back(
//                    std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
//                        multiArcBodiesToPropagate.at( 0 ), multiArcPerBodyPropagationSettings->getInitialStates( ),
//                        arcStartTimesPerBody.at( singleArcBodiesToPropagate.at( i ) ), singleArcBodiesToPropagate.at( i ) ) );

//        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
//        {
//            parameterNamesPerBody.push_back(
//                        std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                            singleArcBodiesToPropagate.at( j ), singleArcInitialState.segment( 6 * j, 6 ), singleArcCentralBodies.at( j ) ) );
//            parameterNamesPerBody.push_back( std::make_shared< EstimatableParameterSettings >(
//                                                 singleArcBodiesToPropagate.at( j ), gravitational_parameter ) );
//        }

//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimatePerBody =
//                createParametersToEstimate< double >( parameterNamesPerBody, bodyMap, hybridArcPerBodyPropagatorSettings );
//        printEstimatableParameterEntries( parametersToEstimatePerBody );

//        HybridArcVariationalEquationsSolver< > perCentralBodyVariationalEquations =
//                HybridArcVariationalEquationsSolver< >(
//                    bodyMap, singleArcIntegratorSettings, multiArcIntegratorSettings,
//                    hybridArcPerBodyPropagatorSettings, parametersToEstimatePerBody, arcStartTimesPerBody.at(
//                        singleArcBodiesToPropagate.at( i ) ), true, false, true/*, dependentVariablesToSave*/ );

//        // Retrieve dependent variables interface.
//        std::shared_ptr< HybridArcDependentVariablesInterface > dependentVariablesInterface =
//                std::dynamic_pointer_cast< HybridArcDependentVariablesInterface >(
//                perCentralBodyVariationalEquations.getDependentVariablesInterface( ) );

//        std::cout << "dependent variables from interface: " << dependentVariablesInterface->getDependentVariables( ( finalTime - initialTime ) / 2.0 ).transpose( ) << "\n\n";

//        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > perBodyMultiArcVariationalSolution =
//                perCentralBodyVariationalEquations.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
//        std::vector< std::map< double, Eigen::VectorXd > > perBodyMultiArcStateSolution =
//                perCentralBodyVariationalEquations.getMultiArcSolver( )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

//        for( unsigned int j = 0; j < perBodyIndicesInFullPropagation.at( singleArcBodiesToPropagate.at( i ) ).size( ); j++ )
//        {
//            for( unsigned int k = 0; k < 2; k++ )
//            {
//                std::map< double, Eigen::MatrixXd > fullMultiArcMatrixHistory = fullMultiArcVariationalSolution.at(
//                            perBodyIndicesInFullPropagation.at( singleArcBodiesToPropagate.at( i ) ).at( j ) ).at( k );
//                std::map< double, Eigen::MatrixXd > perBodyMultiMatrixHistory = perBodyMultiArcVariationalSolution.at( j ).at( k );

//                auto fullIterator = fullMultiArcMatrixHistory.begin( );
//                auto perBodyIterator = perBodyMultiMatrixHistory.begin( );

//                BOOST_CHECK_EQUAL( fullMultiArcMatrixHistory.size( ), perBodyMultiMatrixHistory.size( ) );

//                for( unsigned int i = 0; i < fullMultiArcMatrixHistory.size( ); i++ )
//                {
//                    BOOST_CHECK_CLOSE_FRACTION( fullIterator->first, perBodyIterator->first, std::numeric_limits< double >::epsilon( ) );
//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fullIterator->second, perBodyIterator->second, std::numeric_limits< double >::epsilon( ) );

//                    fullIterator++;
//                    perBodyIterator++;
//                }
//            }
//        }
//    }
}




BOOST_AUTO_TEST_SUITE_END( )

}

}

