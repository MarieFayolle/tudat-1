/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/ppnParameters.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/equivalencePrincipleViolationParameter.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_apparent_distance_partials)

//! Test partial derivatives of apparent distance observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testApparentDistancePartials )
{
    std::cout.precision( 20 );

//    // Define and create ground stations.
//    std::vector< std::pair< std::string, std::string > > groundStations;
//    groundStations.resize( 2 );
//    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
//    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


    // Initial guess central instant.
    double estimatedCentralInstant = 116200.0;

    // Specify initial time
    double initialEphemerisTime = 1.0e7; //estimatedCentralInstant - 3600.0;
    double finalEphemerisTime = 2.0e7; //estimatedCentralInstant + 3600.0;

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Io" );
    bodiesToCreate.push_back( "Europa" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime, finalEphemerisTime );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "ECLIPJ2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                "Earth", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
    bodyMap[ "Earth" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                "Io", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
    bodyMap[ "Io" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                "Europa", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
    bodyMap[ "Europa" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
//    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
//                "Sun", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
//    bodyMap[ "Sun" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );

//    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
//    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
//                "Earth", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
//    bodyMap[ "Earth" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );
//    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
//                "Io", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
//    bodyMap[ "Io" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );
//    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
//                "Europa", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
//    bodyMap[ "Europa" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );

    // Finalize body creation.
//    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair( "Europa" , ""  );



    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
//    parameterPerturbationMultipliers( 2 ) = 10.0;
    // Test partials with constant ephemerides (allows test of position partials)
    {
//        // Create environment
//        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

//        // Set link ends for observation model
//        LinkEnds linkEnds;
//        linkEnds[ transmitter ] = groundStations[ 1 ];
//        linkEnds[ receiver ] = groundStations[ 0 ];

//        // Generate apparent distance model
//        std::vector< std::string > perturbingBodies;
//        perturbingBodies.push_back( "Sun" );
//        std::shared_ptr< ObservationModel< 1 > > apparentDistanceModel =
//                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
//                    linkEnds, std::make_shared< observation_models::ObservationSettings >(
//                        observation_models::apparent_distance, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
//                            perturbingBodies ) ), bodyMap  );

//        // Generate apparent distance model
//        std::vector< std::string > perturbingBodies;
//        perturbingBodies.push_back( "Earth" );
//        std::shared_ptr< ObservationModel< 1 > > apparentDistanceModel =
//                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
//                    linkEnds, std::make_shared< observation_models::ObservationSettings >(
//                        observation_models::apparent_distance/*, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
//                            perturbingBodies )*/ ), bodyMap  );

        // Create light-time correction settings
        std::vector< std::string > lightTimePerturbingBodies = { "Jupiter" };
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                    lightTimePerturbingBodies ) );

        // Create observation settings
        std::shared_ptr< ObservationSettings > observableSettings = std::make_shared< ObservationSettings >
                ( apparent_distance, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
                      ( Eigen::Vector1d( ) << 1.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double > > apparentDistanceModel =
               ObservationModelCreator< 1, double, double >::createObservationModel(
                    linkEnds, observableSettings, bodyMap );


//    // Load spice kernel.
//    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );

//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );
//    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Mars" );

//    // Retrieve body settings
//    std::map< std::string, std::shared_ptr< BodySettings > > defaultBodySettings = getDefaultBodySettings(
//                bodiesToCreate, initialEphemerisTime, finalEphemerisTime );
//    defaultBodySettings[ "Phobos" ] = std::make_shared< BodySettings >( );
//    defaultBodySettings[ "Phobos" ]->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );

//    // Create bodies
//    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

//    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

//    // Define link ends for observations.
//    LinkEnds linkEnds;
//    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
//    linkEnds[ transmitter ] = std::make_pair( "Mars" , ""  );
//    linkEnds[ transmitter2 ] = std::make_pair( "Phobos" , ""  );

//    // Create light-time correction settings
//    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
//    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
//    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
//                                                lightTimePerturbingBodies ) );

//    // Create observation settings
//    std::shared_ptr< ObservationSettings > observableSettings = std::make_shared< ObservationSettings >
//            ( apparent_distance, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
//                  ( Eigen::Vector1d( ) << 1.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

//    // Create observation model.
//    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
//           ObservationModelCreator< 1, double, double >::createObservationModel(
//                linkEnds, observableSettings, bodyMap );
//    std::shared_ptr< ObservationBias< 1 > > observationBias = observationModel->getObservationBiasCalculator( );


        std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatableVectorParameters;
        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedInitialStateParameters;

        estimatedInitialStateParameters.push_back(
                    std::make_shared< InitialTranslationalStateParameter< double > >(
                        "Earth", propagators::getInitialStateOfBody(
                            "Earth", "SSB", bodyMap, 1.1e7 ) ) );
        estimatedInitialStateParameters.push_back(
                    std::make_shared< InitialTranslationalStateParameter< double > >(
                        "Io", propagators::getInitialStateOfBody(
                            "Io", "SSB", bodyMap, 1.1e7 ) ) );
        estimatedInitialStateParameters.push_back(
                    std::make_shared< InitialTranslationalStateParameter< double > >(
                        "Europa", propagators::getInitialStateOfBody(
                            "Europa", "SSB", bodyMap, 1.1e7 ) ) );

        std::shared_ptr< EstimatableParameter< double > > relativisticParameter;
//        relativisticParameter = std::make_shared< estimatable_parameters::PPNParameterGamma >( );
//        relativisticParameter = std::make_shared< estimatable_parameters::EquivalencePrincipleLpiViolationParameter >( );

//        estimatableDoubleParameters.push_back( relativisticParameter );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                std::make_shared< EstimatableParameterSet< double > >( estimatableDoubleParameters,
                        estimatableVectorParameters, estimatedInitialStateParameters );

//        Eigen::Vector1d observationFromReceptionTime = apparentDistanceModel->computeObservations( 1.1e7, receiver );
////        Eigen::Vector1d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData( 1.1e7, receiver, linkEndTimes, linkEndStates );
//        std::vector< Eigen::Vector6d > vectorOfStates;
//        std::vector< double > vectorOfTimes;
//        Eigen::VectorXd currentObservation = apparentDistanceModel->computeObservationsWithLinkEndData(
//                    1.1e7, receiver, vectorOfTimes, vectorOfStates );

//        std::cout << "observation from reception time: " << observationFromReceptionTime << "\n\n";
//        std::cout << "currentObservation: " << currentObservation << "\n\n";
//        for ( int i = 0 ; i < vectorOfTimes.size( ) ; i++ )
//        {
//            std::cout << "times observation: " << vectorOfTimes[ i ] << "\n\n";
//        }

        Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );

//        // Create parameter objects.
//        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
//                createEstimatableParameters( bodyMap, 1.1E7 );

        testObservationPartials( apparentDistanceModel, bodyMap, fullEstimatableParameterSet, linkEnds,
                                 apparent_distance, 1.0E-4, true, true, 1.0, parameterPerturbationMultipliers );
    }


//    // Test partials with real ephemerides (without test of position partials)
//    {
//        std::cout << "Test 1" << std::endl;
//        // Create environment
//        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

//        // Set link ends for observation model
//        LinkEnds linkEnds;
//        linkEnds[ transmitter ] = groundStations[ 1 ];
//        linkEnds[ receiver ] = groundStations[ 0 ];

//        // Generate one-way range model
//        std::shared_ptr< ObservationModel< 2 > > angularPositionModel =
//                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
//                    linkEnds, std::make_shared< observation_models::ObservationSettings >(
//                        observation_models::angular_position ), bodyMap  );

//        // Create parameter objects.
//        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
//                createEstimatableParameters( bodyMap, 1.1E7 );

//        testObservationPartials( angularPositionModel, bodyMap, fullEstimatableParameterSet, linkEnds, angular_position, 1.0E-4,
//        false, true, 1.0, parameterPerturbationMultipliers );

//    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





