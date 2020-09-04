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

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/mutualApproximationObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;


BOOST_AUTO_TEST_SUITE( test_mutual_approximation_model )


BOOST_AUTO_TEST_CASE( testMutualApproximationModel )
{

    // Initial guess central instant.
    double estimatedCentralInstant = 116200.0;

    // Specify initial time
    double initialEphemerisTime = estimatedCentralInstant - 3600.0;
    double finalEphemerisTime = estimatedCentralInstant + 3600.0;

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
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair( "Europa" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< MutualApproximationObservationSettings > observableSettings = std::make_shared< MutualApproximationObservationSettings >
            ( lightTimeCorrectionSettings, 15.0 * 60.0, 30.0, 15.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ), 4, true, false,
              std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector1d( ) << 2.0 ).finished( ), true ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector1d( ) << 5.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ) ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
           ObservationModelCreator< 1, double, double >::createObservationModel( linkEnds, observableSettings, bodyMap );
    std::shared_ptr< ObservationBias< 1 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double receiverObservationTime = estimatedCentralInstant;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector1d observationFromReceptionTime = observationModel->computeObservations( receiverObservationTime, receiver );
    Eigen::Vector1d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimes, linkEndStates );
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 3 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 3 );

    // Manually create and compute light time corrections
    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorFirstTransmitter =
            createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter ], linkEnds[ receiver ] );
    double lightTimeCorrectionFirstTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 2 ), linkEndTimes.at( 0 ), linkEndTimes.at( 2 ) );

    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorSecondTransmitter =
            createLightTimeCorrections(
                lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter2 ], linkEnds[ receiver ] );
    double lightTimeCorrectionSecondTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 1 ), linkEndStates.at( 2 ), linkEndTimes.at( 1 ), linkEndTimes.at( 2 ) );

    // Check equality of computed observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 0 ), bodyMap.at( "Io" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 1 ), bodyMap.at( "Europa" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 2 ), bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 2 ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Check that reception time is kept fixed.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( observationFromReceptionTime[ 0 ] ) - observationBias->getObservationBias(
                                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                                linkEndTimes[ 2 ], std::numeric_limits< double >::epsilon( ) );

    // Manually compute light time
    Eigen::Vector3d positionDifferenceFirstTransmitter = ( linkEndStates[ 0 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceFirstTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionFirstTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 0 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    Eigen::Vector3d positionDifferenceSecondTransmitter = ( linkEndStates[ 1 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceSecondTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionSecondTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 1 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    // Compute apparent distances at required times.
    double observationsStartingTime = estimatedCentralInstant - 15.0 * 60.0;
    std::map< double, double > apparentDistanceObservations;

    // Create apparent distance observation settings
    std::shared_ptr< ObservationSettings > apparentDistanceObservableSettings = std::make_shared< ObservationSettings >
            ( apparent_distance, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector1d( ) << 5.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double > > apparentDistanceObservationModel =
           ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, apparentDistanceObservableSettings, bodyMap );

    // Compute intermediate apparent distance observables.
    double initialisationCounter = - 30.0;
    double counter = initialisationCounter;
    for ( double currentObservationTime = observationsStartingTime ;
          currentObservationTime <= estimatedCentralInstant + 15.0 * 60.0 ; currentObservationTime += 30.0  )
    {
        double currentApparentDistance = apparentDistanceObservationModel->computeObservations( currentObservationTime, receiver )[ 0 ]
                * 3600.0 * 180.0 / mathematical_constants::PI;
        apparentDistanceObservations[ counter ] = currentApparentDistance;
        counter += 1.0;
    }

    // Compute coefficients of the derivative of the fitting polynomial.
    std::vector< double > polynomialPowers = { 0, 1, 2, 3, 4 };
    std::vector< double > polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceObservations, polynomialPowers );

    std::vector< double > derivativePolynomialCoefficients;
    for ( unsigned int i = 1 ; i <= 4 ; i++ )
    {
        derivativePolynomialCoefficients.push_back( i * polynomialCoefficients[ i ] );
    }


    // Compute roots of the derivative of the fitting polynomial.
    double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
    double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
    double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];

    double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
    double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

    double beta = Q * Q * Q + R * R;


    double computedCentralInstant;
    if ( beta < 0 )
    {
        double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
        double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
        double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
        double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;

        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
        if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t2 - initialisationCounter ) * 30.0;
        }
        if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t3 - initialisationCounter ) * 30.0;
        }

    }
    else
    {
        double S = 0.0;
        if ( ( R + sqrt( beta ) ) >= 0 )
        {
            S = std::pow( R + sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            S = - std::pow( std::fabs( R + std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double T = 0.0;
        if ( ( R - sqrt( beta ) ) >= 0 )
        {
            T = std::pow( R - sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            T = - std::pow( std::fabs( R - std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double t1 = - b2 / 3.0 + ( S + T );
        std::cout << "t1: " << t1 << "\n\n";
        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
    }

    BOOST_CHECK_CLOSE_FRACTION(
                computedCentralInstant + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                observationFromReceptionTime( 0 ), std::numeric_limits< double >::epsilon( ) );

    // Compute transmission time from light time.
    double firstTransmitterObservationTime = ( observationFromReceptionTime[ 0 ] - observationBias->getObservationBias(
                std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ) ) - ( linkEndTimes[ 2 ] - linkEndTimes[ 0 ] );
    double secondTransmitterObservationTime = ( observationFromReceptionTime[ 0 ] - observationBias->getObservationBias(
                std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ) ) - ( linkEndTimes[ 2 ] - linkEndTimes[ 1 ] );

    // Compare computed against returned transmission time.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( firstTransmitterObservationTime ), linkEndTimes[ 0 ],
            std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( secondTransmitterObservationTime ), linkEndTimes[ 1 ],
            std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_CASE( testMutualApproximationWithImpactParameterModel )
{

    // Initial guess central instant.
    double estimatedCentralInstant = 116200.0;

    // Specify initial time
    double initialEphemerisTime = estimatedCentralInstant - 3600.0;
    double finalEphemerisTime = estimatedCentralInstant + 3600.0;

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
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair( "Europa" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< MutualApproximationObservationSettings > observableSettings = std::make_shared< MutualApproximationObservationSettings >
            ( lightTimeCorrectionSettings, 15.0 * 60.0, 30.0, 15.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ), 4, true, false,
              std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 2.0, 5.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ) ).finished( ), true ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector1d( ) << 5.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ) ).finished( ), true ),
              true, mutual_approximation_with_impact_parameter );

    // Create observation model.
    std::shared_ptr< ObservationModel< 2, double, double > > observationModel =
           ObservationModelCreator< 2, double, double >::createObservationModel( linkEnds, observableSettings, bodyMap );
    std::shared_ptr< ObservationBias< 2 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double receiverObservationTime = estimatedCentralInstant;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector2d observationFromReceptionTime = observationModel->computeObservations( receiverObservationTime, receiver );
    Eigen::Vector2d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimes, linkEndStates );
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 3 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 3 );

    // Manually create and compute light time corrections
    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorFirstTransmitter =
            createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter ], linkEnds[ receiver ] );
    double lightTimeCorrectionFirstTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 2 ), linkEndTimes.at( 0 ), linkEndTimes.at( 2 ) );

    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorSecondTransmitter =
            createLightTimeCorrections(
                lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter2 ], linkEnds[ receiver ] );
    double lightTimeCorrectionSecondTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 1 ), linkEndStates.at( 2 ), linkEndTimes.at( 1 ), linkEndTimes.at( 2 ) );

    // Check equality of computed observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 0 ), bodyMap.at( "Io" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 1 ), bodyMap.at( "Europa" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 2 ), bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 2 ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Check that reception time is kept fixed.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( observationFromReceptionTime[ 0 ] ) - observationBias->getObservationBias(
                                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                                linkEndTimes[ 2 ], std::numeric_limits< double >::epsilon( ) );

    // Manually compute light time
    Eigen::Vector3d positionDifferenceFirstTransmitter = ( linkEndStates[ 0 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceFirstTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionFirstTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 0 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    Eigen::Vector3d positionDifferenceSecondTransmitter = ( linkEndStates[ 1 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceSecondTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionSecondTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 1 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    // Compute apparent distances at required times.
    double observationsStartingTime = estimatedCentralInstant - 15.0 * 60.0;
    std::map< double, double > apparentDistanceObservations;

    // Create apparent distance observation settings
    std::shared_ptr< ObservationSettings > apparentDistanceObservableSettings = std::make_shared< ObservationSettings >
            ( apparent_distance, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector1d( ) << 5.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double > > apparentDistanceObservationModel =
           ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, apparentDistanceObservableSettings, bodyMap );

    // Compute intermediate apparent distance observables.
    double initialisationCounter = - 30.0;
    double counter = initialisationCounter;
    for ( double currentObservationTime = observationsStartingTime ;
          currentObservationTime <= estimatedCentralInstant + 15.0 * 60.0 ; currentObservationTime += 30.0  )
    {
        double currentApparentDistance = apparentDistanceObservationModel->computeObservations( currentObservationTime, receiver )[ 0 ]
                * 3600.0 * 180.0 / mathematical_constants::PI;
        apparentDistanceObservations[ counter ] = currentApparentDistance;
        counter += 1.0;
    }

    // Compute coefficients of the derivative of the fitting polynomial.
    std::vector< double > polynomialPowers = { 0, 1, 2, 3, 4 };
    std::vector< double > polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceObservations, polynomialPowers );

    std::vector< double > derivativePolynomialCoefficients;
    for ( unsigned int i = 1 ; i <= 4 ; i++ )
    {
        derivativePolynomialCoefficients.push_back( i * polynomialCoefficients[ i ] );
    }


    // Compute roots of the derivative of the fitting polynomial.
    double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
    double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
    double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];

    double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
    double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

    double beta = Q * Q * Q + R * R;


    double computedCentralInstant;
    if ( beta < 0 )
    {
        double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
        double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
        double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
        double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;

        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
        if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t2 - initialisationCounter ) * 30.0;
        }
        if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t3 - initialisationCounter ) * 30.0;
        }

    }
    else
    {
        double S = 0.0;
        if ( ( R + sqrt( beta ) ) >= 0 )
        {
            S = std::pow( R + sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            S = - std::pow( std::fabs( R + std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double T = 0.0;
        if ( ( R - sqrt( beta ) ) >= 0 )
        {
            T = std::pow( R - sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            T = - std::pow( std::fabs( R - std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double t1 = - b2 / 3.0 + ( S + T );
        std::cout << "t1: " << t1 << "\n\n";
        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
    }


    // Compute apparent distance at central instant manually.
    // Compute spherical relative position for first transmitter wrt receiver.
    Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesFirstTransmitter = coordinate_conversions::convertCartesianToSpherical< double >(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ) ).template cast< double >( );

    // Compute spherical relative position for second transmitter wrt receiver.
    Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesSecondTransmitter = coordinate_conversions::convertCartesianToSpherical< double >(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ) ).template cast< double >( );

    // Compute right ascension and declination first transmitter.
    double rightAscensionFirstTransmitter = sphericalRelativeCoordinatesFirstTransmitter.z( );
    double declinationFirstTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesFirstTransmitter.y( );

    // Compute right ascension and declination second transmitter.
    double rightAscensionSecondTransmitter = sphericalRelativeCoordinatesSecondTransmitter.z( );
    double declinationSecondTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesSecondTransmitter.y( );

    // Compute partials of right ascension w.r.t. time.
    double timeDerivativeRightAscensionFirstTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 0 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );
    double timeDerivativeRightAscensionSecondTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 1 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );

    // Compute partials of declination w.r.t. time.
    double timeDerivativeDeclinationFirstTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 0 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );
    double timeDerivativeDeclinationSecondTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 1 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );

    // Compute average declination.
    double averageDeclination = ( declinationFirstTransmitter + declinationSecondTransmitter ) / 2.0;

    Eigen::Vector2d relativePositionInReceiverFrame =
            ( Eigen::Vector2d( ) << ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter ) * std::cos( averageDeclination ),
              declinationSecondTransmitter - declinationFirstTransmitter ).finished( );

    Eigen::Vector2d relativeVelocityInReceiverFrame =
            ( Eigen::Vector2d( ) << ( timeDerivativeRightAscensionSecondTransmitter - timeDerivativeRightAscensionFirstTransmitter )
            * std::cos( averageDeclination )
            - ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter ) / 2.0
            * std::sin( averageDeclination ) * ( timeDerivativeDeclinationFirstTransmitter + timeDerivativeDeclinationSecondTransmitter ),
              timeDerivativeDeclinationSecondTransmitter - timeDerivativeDeclinationFirstTransmitter ).finished( );

    double apparentDistance = std::sqrt( relativePositionInReceiverFrame[ 0 ] * relativePositionInReceiverFrame[ 0 ]
            + relativePositionInReceiverFrame[ 1 ] * relativePositionInReceiverFrame[ 1 ] );

    // Compute modified observable.
    double apparentDistanceAtCentralInstant = 1.0 / ( apparentDistance )
            * ( relativePositionInReceiverFrame[ 0 ] * relativeVelocityInReceiverFrame[ 0 ]
            + relativePositionInReceiverFrame[ 1 ] * relativeVelocityInReceiverFrame[ 1 ] );

    BOOST_CHECK_CLOSE_FRACTION( computedCentralInstant + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                observationFromReceptionTime( 0 ), std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( apparentDistance + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).y( ),
                observationFromReceptionTime( 1 ), std::numeric_limits< double >::epsilon( ) );

    // Compute transmission time from light time.
    double firstTransmitterObservationTime = ( observationFromReceptionTime[ 0 ] - observationBias->getObservationBias(
                std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ) ) - ( linkEndTimes[ 2 ] - linkEndTimes[ 0 ] );
    double secondTransmitterObservationTime = ( observationFromReceptionTime[ 0 ] - observationBias->getObservationBias(
                std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ) ) - ( linkEndTimes[ 2 ] - linkEndTimes[ 1 ] );

    // Compare computed against returned transmission time.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( firstTransmitterObservationTime ), linkEndTimes[ 0 ],
            std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( secondTransmitterObservationTime ), linkEndTimes[ 1 ],
            std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_CASE( testModifiedMutualApproximationModel )
{

    // Initial guess central instant.
    double estimatedCentralInstant = 116200.0;

    // Specify initial time
    double initialEphemerisTime = estimatedCentralInstant - 3600.0;
    double finalEphemerisTime = estimatedCentralInstant + 3600.0;

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
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair( "Europa" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< MutualApproximationObservationSettings > observableSettings = std::make_shared< MutualApproximationObservationSettings >
            ( lightTimeCorrectionSettings, 15.0 * 60.0, 30.0, 15.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ), 4, true, false,
              std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector1d( ) << 5.0e-13 ).finished( ), true ),
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector1d( ) << 5.0 * mathematical_constants::PI / ( 3600.0 * 180.0 ) ).finished( ), true ),
              false );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
           ObservationModelCreator< 1, double, double >::createObservationModel( linkEnds, observableSettings, bodyMap );
    std::shared_ptr< ObservationBias< 1 > > observationBias = observationModel->getObservationBiasCalculator( );


    // Compute observation separately with two functions.
    double receiverObservationTime = estimatedCentralInstant;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector1d observationFromReceptionTime = observationModel->computeObservations( receiverObservationTime, receiver );
    Eigen::Vector1d observationFromReceptionTime2 = observationModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimes, linkEndStates );
    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 3 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 3 );

    // Manually create and compute light time corrections
    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorFirstTransmitter =
            createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter ], linkEnds[ receiver ] );
    double lightTimeCorrectionFirstTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 2 ), linkEndTimes.at( 0 ), linkEndTimes.at( 2 ) );

    std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculatorSecondTransmitter =
            createLightTimeCorrections(
                lightTimeCorrectionSettings.at( 0 ), bodyMap, linkEnds[ transmitter2 ], linkEnds[ receiver ] );
    double lightTimeCorrectionSecondTransmitter = lightTimeCorrectionCalculatorFirstTransmitter->calculateLightTimeCorrection(
                linkEndStates.at( 1 ), linkEndStates.at( 2 ), linkEndTimes.at( 1 ), linkEndTimes.at( 2 ) );

    // Check equality of computed observations.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    // Check consistency of link end states and time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 0 ), bodyMap.at( "Io" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 1 ), bodyMap.at( "Europa" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStates.at( 2 ), bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 2 ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Manually compute light time
    Eigen::Vector3d positionDifferenceFirstTransmitter = ( linkEndStates[ 0 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceFirstTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionFirstTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 0 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    Eigen::Vector3d positionDifferenceSecondTransmitter = ( linkEndStates[ 1 ] - linkEndStates[ 2 ] ).segment( 0, 3 );
    BOOST_CHECK_CLOSE_FRACTION(
                positionDifferenceSecondTransmitter.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrectionSecondTransmitter,
                linkEndTimes[ 2 ] - linkEndTimes[ 1 ],
                std::numeric_limits< double >::epsilon( ) * 1000.0 );
                        // Poor tolerance due to rounding errors when subtracting times

    // Compute apparent distances at required times.
    double observationsStartingTime = estimatedCentralInstant - 15.0 * 60.0;
    std::map< double, double > apparentDistanceObservations;

    // Create apparent distance observation settings
    std::shared_ptr< ObservationSettings > apparentDistanceObservableSettings = std::make_shared< ObservationSettings >
            ( apparent_distance, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
                  ( Eigen::Vector1d( ) << 5.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double > > apparentDistanceObservationModel =
           ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, apparentDistanceObservableSettings, bodyMap );

    // Compute intermediate apparent distance observables.
    double initialisationCounter = - 30.0;
    double counter = initialisationCounter;
    for ( double currentObservationTime = observationsStartingTime ;
          currentObservationTime <= estimatedCentralInstant + 15.0 * 60.0 ; currentObservationTime += 30.0  )
    {
        double currentApparentDistance = apparentDistanceObservationModel->computeObservations( currentObservationTime, receiver )[ 0 ]
                * 3600.0 * 180.0 / mathematical_constants::PI;
        apparentDistanceObservations[ counter ] = currentApparentDistance;
        counter += 1.0;
    }

    // Compute coefficients of the derivative of the fitting polynomial.
    std::vector< double > polynomialPowers = { 0, 1, 2, 3, 4 };
    std::vector< double > polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceObservations, polynomialPowers );

    std::vector< double > derivativePolynomialCoefficients;
    for ( unsigned int i = 1 ; i <= 4 ; i++ )
    {
        derivativePolynomialCoefficients.push_back( i * polynomialCoefficients[ i ] );
    }


    // Compute roots of the derivative of the fitting polynomial.
    double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
    double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
    double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];

    double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
    double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

    double beta = Q * Q * Q + R * R;


    double computedCentralInstant;
    if ( beta < 0 )
    {
        double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
        double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
        double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
        double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;

        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
        if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t2 - initialisationCounter ) * 30.0;
        }
        if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t3 - initialisationCounter ) * 30.0;
        }

    }
    else
    {
        double S = 0.0;
        if ( ( R + sqrt( beta ) ) >= 0 )
        {
            S = std::pow( R + sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            S = - std::pow( std::fabs( R + std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double T = 0.0;
        if ( ( R - sqrt( beta ) ) >= 0 )
        {
            T = std::pow( R - sqrt( beta ), 1.0 / 3.0 );
        }
        else
        {
            T = - std::pow( std::fabs( R - std::sqrt( beta ) ), 1.0 / 3.0 );
        }

        double t1 = - b2 / 3.0 + ( S + T );
        std::cout << "t1: " << t1 << "\n\n";
        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            computedCentralInstant = observationsStartingTime + ( t1 - initialisationCounter ) * 30.0;
        }
    }


    // Compute apparent distance at central instant manually.
    // Compute spherical relative position for first transmitter wrt receiver.
    Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesFirstTransmitter = coordinate_conversions::convertCartesianToSpherical< double >(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ) ).template cast< double >( );

    // Compute spherical relative position for second transmitter wrt receiver.
    Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesSecondTransmitter = coordinate_conversions::convertCartesianToSpherical< double >(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ) ).template cast< double >( );

    // Compute right ascension and declination first transmitter.
    double rightAscensionFirstTransmitter = sphericalRelativeCoordinatesFirstTransmitter.z( );
    double declinationFirstTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesFirstTransmitter.y( );

    // Compute right ascension and declination second transmitter.
    double rightAscensionSecondTransmitter = sphericalRelativeCoordinatesSecondTransmitter.z( );
    double declinationSecondTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesSecondTransmitter.y( );

    // Compute partials of right ascension w.r.t. time.
    double timeDerivativeRightAscensionFirstTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 0 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );
    double timeDerivativeRightAscensionSecondTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 1 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );

    // Compute partials of declination w.r.t. time.
    double timeDerivativeDeclinationFirstTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                linkEndStates[ 0 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 0 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );
    double timeDerivativeDeclinationSecondTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 2 ].segment( 0, 3 ),
                linkEndStates[ 1 ].segment( 3, 3 ) - linkEndStates[ 2 ].segment( 3, 3 ) );

    // Compute average declination.
    double averageDeclination = ( declinationFirstTransmitter + declinationSecondTransmitter ) / 2.0;

    Eigen::Vector2d relativePositionInReceiverFrame =
            ( Eigen::Vector2d( ) << ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter ) * std::cos( averageDeclination ),
              declinationSecondTransmitter - declinationFirstTransmitter ).finished( );

    Eigen::Vector2d relativeVelocityInReceiverFrame =
            ( Eigen::Vector2d( ) << ( timeDerivativeRightAscensionSecondTransmitter - timeDerivativeRightAscensionFirstTransmitter )
            * std::cos( averageDeclination )
            - ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter ) / 2.0
            * std::sin( averageDeclination ) * ( timeDerivativeDeclinationFirstTransmitter + timeDerivativeDeclinationSecondTransmitter ),
              timeDerivativeDeclinationSecondTransmitter - timeDerivativeDeclinationFirstTransmitter ).finished( );

    double apparentDistance = std::sqrt( relativePositionInReceiverFrame[ 0 ] * relativePositionInReceiverFrame[ 0 ]
            + relativePositionInReceiverFrame[ 1 ] * relativePositionInReceiverFrame[ 1 ] );

    // Compute modified observable.
    double modifiedMutualApproximationObservable = 1.0 / ( apparentDistance )
            * ( relativePositionInReceiverFrame[ 0 ] * relativeVelocityInReceiverFrame[ 0 ]
            + relativePositionInReceiverFrame[ 1 ] * relativeVelocityInReceiverFrame[ 1 ] );

    BOOST_CHECK_CLOSE_FRACTION( modifiedMutualApproximationObservable + observationBias->getObservationBias(
                    std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                observationFromReceptionTime( 0 ), std::numeric_limits< double >::epsilon( ) );

    // Check that reception time is kept fixed.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( computedCentralInstant ), linkEndTimes[ 2 ], std::numeric_limits< double >::epsilon( ) );

    // Compute transmission time from light time.
    double firstTransmitterObservationTime = computedCentralInstant - ( linkEndTimes[ 2 ] - linkEndTimes[ 0 ] );
    double secondTransmitterObservationTime = computedCentralInstant - ( linkEndTimes[ 2 ] - linkEndTimes[ 1 ] );

    // Compare computed against returned transmission time.
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( firstTransmitterObservationTime ), linkEndTimes[ 0 ],
            std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( secondTransmitterObservationTime ), linkEndTimes[ 1 ],
            std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_SUITE_END( )

}

}

