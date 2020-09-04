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


double computePartialRightAscensionWrtTime( Eigen::Vector3d cartesianPosition,
                                            Eigen::Vector3d cartesianVelocity )
{
    double rx = cartesianPosition[ 0 ];
    double ry = cartesianPosition[ 1 ];
    double derivRx = cartesianVelocity[ 0 ];
    double derivRy = cartesianVelocity[ 1 ];

    double partial = 2.0 * ( std::sqrt( rx * rx + ry * ry ) + rx )
            / ( ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx ) + ry * ry ) * 1.0 / std::sqrt( rx * rx + ry * ry )
            * ( rx * derivRy - ry * derivRx );

    return partial;
}


double computePartialDeclinationWrtTime( Eigen::Vector3d cartesianPosition,
                                         Eigen::Vector3d cartesianVelocity )
{
    double rx = cartesianPosition[ 0 ];
    double ry = cartesianPosition[ 1 ];
    double rz = cartesianPosition[ 2 ];
    double derivRx = cartesianVelocity[ 0 ];
    double derivRy = cartesianVelocity[ 1 ];
    double derivRz = cartesianVelocity[ 2 ];

    double squaredDistance = rx * rx + ry * ry + rz * rz;

    double partial = ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * derivRz
            - ( rz / ( std::sqrt( rx * rx + ry * ry ) * squaredDistance ) ) * ( rx * derivRx + ry * derivRy + rz * derivRz );

    return partial;
}

double computeSecondPartialXwrtTime( Eigen::Vector3d cartesianPositionFirstObject,
                                     Eigen::Vector3d cartesianPositionSecondObject,
                                     Eigen::Vector3d cartesianVelocityFirstObject,
                                     Eigen::Vector3d cartesianVelocitySecondObject,
                                     Eigen::Vector3d cartesianAccelerationFirstObject,
                                     Eigen::Vector3d cartesianAccelerationSecondObject,
                                     double rightAscensionFirstObject, double declinationFirstObject,
                                     double rightAscensionSecondObject, double declinationSecondObject )
{
    double averageDeclination = ( declinationFirstObject + declinationSecondObject ) / 2.0;

    double partialRightAscensionFirstObject = computePartialRightAscensionWrtTime( cartesianPositionFirstObject,
                                                                                   cartesianVelocityFirstObject );
    double partialRightAscensionSecondObject = computePartialRightAscensionWrtTime( cartesianPositionSecondObject,
                                                                                    cartesianVelocitySecondObject );

    double partialDeclinationFirstObject = computePartialDeclinationWrtTime( cartesianPositionFirstObject,
                                                                             cartesianVelocityFirstObject );
    double partialDeclinationSecondObject = computePartialDeclinationWrtTime( cartesianPositionSecondObject,
                                                                              cartesianVelocitySecondObject );

    double secondPartialRightAscensionFirstObject = computeSecondPartialRightAscensionWrtTime( cartesianPositionFirstObject,
                                                                                               cartesianVelocityFirstObject,
                                                                                               cartesianAccelerationFirstObject );
    double secondPartialRightAscensionSecondObject = computeSecondPartialRightAscensionWrtTime( cartesianPositionSecondObject,
                                                                                                cartesianVelocitySecondObject,
                                                                                                cartesianAccelerationSecondObject );

    double secondPartialDeclinationFirstObject = computeSecondPartialDeclinationWrtTime( cartesianPositionFirstObject,
                                                                                         cartesianVelocityFirstObject,
                                                                                         cartesianAccelerationFirstObject );
    double secondPartialDeclinationSecondObject = computeSecondPartialDeclinationWrtTime( cartesianPositionSecondObject,
                                                                                          cartesianVelocitySecondObject,
                                                                                          cartesianAccelerationSecondObject );

    double partial = ( secondPartialRightAscensionSecondObject - secondPartialRightAscensionFirstObject ) * std::cos( averageDeclination )
            - ( partialRightAscensionSecondObject - partialRightAscensionFirstObject ) * std::sin( averageDeclination )
            * ( partialDeclinationFirstObject + partialDeclinationSecondObject )
            - ( ( rightAscensionSecondObject - rightAscensionFirstObject ) / 4.0 ) * ( partialDeclinationFirstObject + partialDeclinationSecondObject )
            * ( partialDeclinationFirstObject + partialDeclinationSecondObject ) * std::cos( averageDeclination )
            - ( ( rightAscensionSecondObject - rightAscensionFirstObject ) / 2.0 ) * std::sin( averageDeclination ) *
            ( secondPartialDeclinationFirstObject + secondPartialDeclinationSecondObject );

    return partial;
}


double computeSecondPartialYwrtTime( Eigen::Vector3d cartesianPositionFirstObject,
                                     Eigen::Vector3d cartesianPositionSecondObject,
                                     Eigen::Vector3d cartesianVelocityFirstObject,
                                     Eigen::Vector3d cartesianVelocitySecondObject,
                                     Eigen::Vector3d cartesianAccelerationFirstObject,
                                     Eigen::Vector3d cartesianAccelerationSecondObject )
{
    double secondPartialDeclinationFirstObject = computeSecondPartialDeclinationWrtTime( cartesianPositionFirstObject,
                                                                                         cartesianVelocityFirstObject,
                                                                                         cartesianAccelerationFirstObject );
    double secondPartialDeclinationSecondObject = computeSecondPartialDeclinationWrtTime( cartesianPositionSecondObject,
                                                                                          cartesianVelocitySecondObject,
                                                                                          cartesianAccelerationSecondObject );

    double partial = secondPartialDeclinationSecondObject - secondPartialDeclinationFirstObject;

    return partial;
}


std::pair< double, double > computeRightAscensionDeclination( Eigen::Vector3d relativePosition )
{
    Eigen::Matrix< double, 3, 1 > sphericalCoordinates = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                relativePosition ).template cast< double >( );
    double rightAscension = sphericalCoordinates.z( );
    double declination = mathematical_constants::PI / 2.0 - sphericalCoordinates.y( );

    return std::make_pair( rightAscension, declination );
}

double computeCentralInstantFromDependentVariables(
        double simulationStartEpoch,
        double frequencyApparentDistanceObservations,
        std::map< double, Eigen::VectorXd > dependentVariablesHistory )
{
    Eigen::Matrix< double, Eigen::Dynamic, 1 > apparentDistances;
    apparentDistances.resize( dependentVariablesHistory.size( ) );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > apparentDistanceTimes;
    apparentDistanceTimes.resize( dependentVariablesHistory.size( ) );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > realApparentDistanceTimes;
    realApparentDistanceTimes.resize( dependentVariablesHistory.size( ) );

    double initialisationCounter = - ( int ) ( dependentVariablesHistory.size( ) / 2.0 );
//    std::cout << "TEST: " << - ( int ) ( dependentVariablesHistory.size( ) / 2.0 ) << "\n\n";

    double observationStartingTime = simulationStartEpoch;

    int counter = 0;
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ; itr != dependentVariablesHistory.end( ) ; itr++ )
    {
        std::pair< double, double > currentRightAscensionDeclinationIo = computeRightAscensionDeclination( itr->second.segment( 9, 3 ) );
        std::pair< double, double > currentRightAscensionDeclinationEuropa = computeRightAscensionDeclination( itr->second.segment( 15, 3 ) );
        double currentAverageDeclination = ( currentRightAscensionDeclinationIo.second + currentRightAscensionDeclinationEuropa.second ) / 2.0;
        apparentDistanceTimes( counter, 0 ) = counter + initialisationCounter; // itr->first;
        realApparentDistanceTimes( counter, 0 ) = itr->first;
        apparentDistances( counter, 0 ) = std::sqrt(
                    ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
                    * ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
                    + ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second )
                    * ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second ) ) * 3600.0 * 180.0 / mathematical_constants::PI;
        counter += 1;
//        std::cout << "apparent distance Io-Europa (arcseconds): " << std::sqrt(
//                         ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                         * ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                         + ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second )
//                         * ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second ) ) * 3600.0 * 180.0 / mathematical_constants::PI << "at time = " << itr->first << " seconds." "\n\n";
    }

//    std::cout << "apparentDistanceTimes: " << realApparentDistanceTimes.transpose( ) << "\n\n";
//    std::cout << "apparentDistances: " << apparentDistances.transpose( ) << "\n\n";

    std::vector< double > polynomialPowers = { 0, 1, 2, 3, 4 }; //, 4 };
    Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceTimes, apparentDistances, polynomialPowers );
    std::cout << "polynomial coefficients: " << polynomialCoefficients.transpose( ) << "\n\n";
    Eigen::VectorXd derivativePolynomialCoefficients; derivativePolynomialCoefficients.resize( polynomialCoefficients.size( ) - 1 );
    for ( unsigned int i = 1 ; i <= 4 ; i++ )
    {
        derivativePolynomialCoefficients[ i - 1 ] = i * polynomialCoefficients[ i ];
    }


    double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
    double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
    double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];

    double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
    double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

    double beta = Q * Q * Q + R * R;



    double estimatedCentralInstant;
    if ( beta < 0 )
    {
        double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
        double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
        double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
        double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
//        std::cout << "t1: " << t1 << "\n\n";
//        std::cout << "t2: " << t2 << "\n\n";
//        std::cout << "t3: " << t3 << "\n\n";
//        std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
        if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
        {
            estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations;
//            std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";
        }
        if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
        {
            estimatedCentralInstant = observationStartingTime + ( t2 - initialisationCounter ) * frequencyApparentDistanceObservations;
//            std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";
        }
        if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
        {
            estimatedCentralInstant = observationStartingTime + ( t3 - initialisationCounter ) * frequencyApparentDistanceObservations;
//            std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";
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
            estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations;
        }
    }

    return estimatedCentralInstant;
}


Eigen::Vector3d calculatePointMassGravityAcceleration(
        const Eigen::Vector3d relativePosition,
        double gravitationalParameter )
{
    double relativePositionNorm = relativePosition.norm( );
    Eigen::Vector3d accelerationVector = gravitationalParameter * relativePosition
            / ( relativePositionNorm * relativePositionNorm * relativePositionNorm );

    return accelerationVector;
}

// Function that checks whether the apparent relative distance between two objects in a list, as seen from a common observer, meets a certain requirement.
bool isApparentDistanceBelowThreshold( const double time, const double threshold, const std::string observer, const std::vector< std::string > listOfObjects,
                                       tudat::simulation_setup::NamedBodyMap& bodyMap, const double timeLimit, bool useThresholdAsLowerLimit )
{
    bool apparentDistanceBelowThreshold = false;

    unsigned int numberOfObjects = listOfObjects.size( );

    Eigen::Vector6d cartesianStateObserver = bodyMap.at( observer )->getState( );

    for ( unsigned int firstObjectIndex = 0 ; firstObjectIndex < numberOfObjects - 1 ; firstObjectIndex++ )
    {
        for ( unsigned int secondObjectIndex = firstObjectIndex + 1 ; secondObjectIndex < numberOfObjects ; secondObjectIndex++ )
        {
            std::string firstObject = listOfObjects[ firstObjectIndex ];
            std::string secondObject = listOfObjects[ secondObjectIndex ];

            Eigen::Vector6d cartesianStateFirstObject = bodyMap.at( firstObject )->getState( );
            Eigen::Vector6d cartesianStateSecondObject = bodyMap.at( secondObject )->getState( );

            Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesFirstObject = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                        cartesianStateFirstObject.segment( 0, 3 ) - cartesianStateObserver.segment( 0, 3 ) ).template cast< double >( );
            Eigen::Matrix< double, 3, 1 > sphericalRelativeCoordinatesSecondObject = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                        cartesianStateSecondObject.segment( 0, 3 ) - cartesianStateObserver.segment( 0, 3 ) ).template cast< double >( );

            double rightAscensionFirstObject = sphericalRelativeCoordinatesFirstObject.z( );
            double declinationFirstObject = tudat::mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesFirstObject.y( );

            double rightAscensionSecondObject = sphericalRelativeCoordinatesSecondObject.z( );
            double declinationSecondObject = tudat::mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesSecondObject.y( );

            double differenceRightAscension = rightAscensionSecondObject - rightAscensionFirstObject;
            double differenceDeclination = declinationSecondObject - declinationFirstObject;

            double apparentDistance = std::sqrt( ( differenceRightAscension * std::cos( declinationFirstObject ) ) *
                                                 ( differenceRightAscension * std::cos( declinationFirstObject ) )
                                                 + differenceDeclination * differenceDeclination );

            if ( useThresholdAsLowerLimit )
            {
                if ( apparentDistance <= threshold )
                {
                    std::cout << "apparent distance between " << firstObject << " and " << secondObject << ": "
                              << apparentDistance * 180.0 / tudat::mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
                    std::cout << "threshold: " << threshold * 180.0 / tudat::mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
                    apparentDistanceBelowThreshold = true;
                }
            }
            else
            {
                if ( apparentDistance > threshold )
                {
                    std::cout << "apparent distance between " << firstObject << " and " << secondObject << ": "
                              << apparentDistance * 180.0 / tudat::mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
                    std::cout << "threshold: " << threshold * 180.0 / tudat::mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
                    apparentDistanceBelowThreshold = true;
                }
            }

        }
    }

    if ( time >= timeLimit )
    {
        apparentDistanceBelowThreshold = true;
    }

    return apparentDistanceBelowThreshold;

}


BOOST_AUTO_TEST_SUITE( test_mutual_approximation_partials)

//! Test partial derivatives of mutual approximation observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testIntermediateVariablesOfMutualApproximationPartials )
{

    std::cout.precision( 16 );

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 1341680.0 - 0.5 * 3600.0; // 0.0;
    const double simulationEndEpoch = 1341680.0 + 0.5 * 3600.0; // 1.0 * tudat::physical_constants::JULIAN_YEAR; // 1.0 / 24.0 * tudat::physical_constants::JULIAN_DAY; // 2.0 * tudat::physical_constants::JULIAN_YEAR;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Io" );
    bodiesToCreate.push_back( "Europa" );
//    bodiesToCreate.push_back( "Ganymede" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsEarth;
    currentAccelerationsEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsJupiter;
    currentAccelerationsJupiter[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsJupiter[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsJupiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsIo;
    currentAccelerationsIo[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsEuropa;
    currentAccelerationsEuropa[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );


    accelerationMap[ "Earth" ] = currentAccelerationsEarth;
    accelerationMap[ "Jupiter" ] = currentAccelerationsJupiter;
    accelerationMap[ "Io" ] = currentAccelerationsIo;
    accelerationMap[ "Europa" ] = currentAccelerationsEuropa;

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Earth" );
    bodiesToPropagate.push_back( "Jupiter" );
    bodiesToPropagate.push_back( "Io" );
    bodiesToPropagate.push_back( "Europa" );

    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "SSB" );
    centralBodies.push_back( "SSB" );
    centralBodies.push_back( "Jupiter" );
    centralBodies.push_back( "Jupiter" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = propagators::getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, simulationStartEpoch );
//    systemInitialState[ 7 ] *= 0.999;


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;

//    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList;
//    terminationSettingsList.push_back( std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

//    double limitValueImpactParameter = 10.0 / 3600.0 * mathematical_constants::PI / 180.0; //5.0 * physical_constants::ASTRONOMICAL_UNIT * 70.0 / 3600.0 * mathematical_constants::PI / 180.0;
//    std::cout << "limit value impact parameter: " << limitValueImpactParameter << "\n\n";

    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Io" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Europa" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Io", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Io", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Europa", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Europa", "Earth" ) );

    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Io", "Jupiter" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Earth", "Jupiter" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Europa", "Jupiter" ) );

    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Europa", "Jupiter" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Earth", "Jupiter" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Io", "Jupiter" ) );

    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Earth", "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Io", "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Europa", "Sun" ) );

    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Jupiter" ) );

    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Earth", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Earth", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Jupiter", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Jupiter", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Io", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Io", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_position_dependent_variable, "Europa", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::relative_velocity_dependent_variable, "Europa", "Sun" ) );


    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Define propagator settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                std::make_shared< propagators::PropagationTimeTerminationSettings > ( simulationEndEpoch, true ), propagators::cowell, dependentVariablesToSave );

    // Define numerical integrator settings.
    std::shared_ptr< numerical_integrators:: IntegratorSettings< > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< > >( numerical_integrators::rungeKutta4, simulationStartEpoch, 60.0 /*3600.0*/ ); //0.5 ); //3600.0 / 10.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBITS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create simulation object and propagate dynamics.
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Earth", systemInitialState.segment( 0, 6 ), "SSB" ) );
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Jupiter", systemInitialState.segment( 6, 6 ), "SSB" ) );
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Io", systemInitialState.segment( 12, 6 ), "Jupiter" ) );
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Europa", systemInitialState.segment( 18, 6 ), "Jupiter" ) );


    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    // Create simulation object and propagate dynamics.
    propagators::SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
                dependentVariablesToSave );

    std::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator = variationalEquationsSimulator.getDynamicsSimulator( );

    std::map< double, Eigen::MatrixXd > stateTransitionHistory = variationalEquationsSimulator.getNumericalVariationalEquationsSolution( )[ 0 ];
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( stateTransitionHistory ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( stateTransitionHistory ), 4 );


    // Create state history interpolator.
    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > stateHistoryInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( stateHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( stateHistory ), 4 );

    // Retrieve dependent variables history.
    std::map< double, Eigen::VectorXd > dependentVariablesHistory = dynamicsSimulator->getDependentVariableHistory( );
    // Create dependent variables interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory ), 4 );

    std::shared_ptr< propagators::SingleArcDependentVariablesInterface > dependentVariablesInterface = std::make_shared< propagators::SingleArcDependentVariablesInterface >(
                dependentVariablesInterpolator, dependentVariablesToSave );

    std::shared_ptr< observation_partials::MutualApproximationScaling > mutualApproximationScaling =
            std::make_shared< observation_partials::MutualApproximationScaling >( dependentVariablesInterface );

    std::shared_ptr< observation_partials::MutualApproximationScaling > modifiedMutualApproximationScaling =
            std::make_shared< observation_partials::MutualApproximationScaling >( dependentVariablesInterface );

    double estimatedCentralInstant = computeCentralInstantFromDependentVariables( simulationStartEpoch, 60.0,
                                                                                  dependentVariablesHistory );
    std::cout << "estimatedCentralInstant: " << estimatedCentralInstant << "\n\n";
//    std::cout << "difference w.r.t. matlab: " << estimatedCentralInstant - 1.341684875158209e+06 << "\n\n";



//    double timeVariation = 0.5; //3600.0 / 10.0;

//    std::map< double, Eigen::VectorXd > dependentVariables = dynamicsSimulator->getDependentVariableHistory( );
//    std::map< double, Eigen::VectorXd >::iterator itr = stateHistory.end( );
//    itr--;
//    itr--;
//    double finalPropagationEpoch = itr->first;
//    Eigen::VectorXd finalPropagatedState = itr->second;
//    Eigen::VectorXd finalStateIo = finalPropagatedState.segment( 12, 6 );
//    Eigen::VectorXd finalStateEuropa = finalPropagatedState.segment( 18, 6 );
//    Eigen::VectorXd accelerationIo = dependentVariables[ itr->first ].segment( 0, 3 );
//    Eigen::VectorXd accelerationEuropa = dependentVariables[ itr->first ].segment( 3, 3 );
//    Eigen::VectorXd accelerationEarth = dependentVariables[ itr->first ].segment( 6, 3 );
//    Eigen::VectorXd finalRelativeStateIoEarth = dependentVariables[ itr->first ].segment( 9, 6 );
//    Eigen::VectorXd finalRelativeStateEuropaEarth = dependentVariables[ itr->first ].segment( 15, 6 );
//    Eigen::Vector3d finalRelativeAccelerationIoEarth = ( accelerationIo - accelerationEarth ).segment( 0, 3 );

//    for ( int i = 0 ; i < 1 ; i++ ){
//    itr--;}
//    double preFinalPropagationEpoch = itr->first;
//    Eigen::VectorXd preFinalPropagatedState = itr->second;
//    Eigen::VectorXd preFinalStateIo = preFinalPropagatedState.segment( 12, 6 );
//    Eigen::VectorXd preFinalStateEuropa = preFinalPropagatedState.segment( 18, 6 );
//    Eigen::VectorXd preFinalAccelerationIo = dependentVariables[ itr->first ].segment( 0, 3 );
//    Eigen::VectorXd preFinalAccelerationEuropa = dependentVariables[ itr->first ].segment( 3, 3 );
//    Eigen::VectorXd preFinalAccelerationEarth = dependentVariables[ itr->first ].segment( 6, 3 );
//    Eigen::VectorXd preFinalRelativeStateIoEarth = dependentVariables[ itr->first ].segment( 9, 6 );
//    Eigen::VectorXd preFinalRelativeStateEuropaEarth = dependentVariables[ itr->first ].segment( 15, 6 );
//    Eigen::Vector3d preFinalRelativeAccelerationIoEarth = ( preFinalAccelerationIo - preFinalAccelerationEarth ).segment( 0, 3 );

//    for ( int i = 0 ; i < 1 ; i++ ){
//    itr--;}
//    double prepreFinalPropagationEpoch = itr->first;
//    Eigen::VectorXd prepreFinalPropagatedState = itr->second;
//    Eigen::VectorXd prepreFinalStateIo = prepreFinalPropagatedState.segment( 12, 6 );
//    Eigen::VectorXd prepreFinalStateEuropa = prepreFinalPropagatedState.segment( 18, 6 );
//    Eigen::VectorXd prepreFinalAccelerationIo = dependentVariables[ itr->first ].segment( 0, 3 );
//    Eigen::VectorXd prepreFinalAccelerationEuropa = dependentVariables[ itr->first ].segment( 3, 3 );
//    Eigen::VectorXd prepreFinalAccelerationEarth = dependentVariables[ itr->first ].segment( 6, 3 );
//    Eigen::VectorXd prepreFinalRelativeStateIoEarth = dependentVariables[ itr->first ].segment( 9, 6 );
//    Eigen::VectorXd prepreFinalRelativeStateEuropaEarth = dependentVariables[ itr->first ].segment( 15, 6 );
//    Eigen::Vector3d prepreFinalRelativeAccelerationIoEarth = ( prepreFinalAccelerationIo - prepreFinalAccelerationEarth ).segment( 0, 3 );

////    std::cout << "position vector Io: " << ( finalStateIo - preFinalStateIo ).segment( 0, 3 ).transpose( ) << "\n\n";
////    std::cout << "approximate position Io: " << timeVariation * preFinalStateIo.segment( 3, 3 ).transpose( ) << "\n\n";

////    std::cout << "velocity vector Io: " << ( finalStateIo - preFinalStateIo ).segment( 3, 3 ).transpose( ) << "\n\n";
////    std::cout << "approximate velocity Io: " << timeVariation * preFinalAccelerationIo.transpose( ) << "\n\n";

//    timeVariation = finalPropagationEpoch - preFinalPropagationEpoch;
//    double preTimeVariation = preFinalPropagationEpoch - prepreFinalPropagationEpoch;


//    // Update mutual approximation scaling.
//    std::vector< Eigen::Vector6d > linkEndStates;
//    linkEndStates.push_back( preFinalRelativeStateIoEarth );
//    linkEndStates.push_back( preFinalRelativeStateEuropaEarth );
//    linkEndStates.push_back( Eigen::Vector6d::Zero( ) );

//    std::vector< double > times;
//    times.push_back( preFinalPropagationEpoch );
//    times.push_back( preFinalPropagationEpoch );
//    times.push_back( preFinalPropagationEpoch );

//    observation_models::LinkEnds linkEnds;
//    linkEnds[ observation_models::receiver ] = std::make_pair( "Earth" , ""  );
//    linkEnds[ observation_models::transmitter ] = std::make_pair( "Io" , ""  );
//    linkEnds[ observation_models::transmitter2 ] = std::make_pair( "Europa" , ""  );


//    Eigen::VectorXd currentObservation;
//    mutualApproximationScaling->update( linkEndStates, times, observation_models::receiver, linkEnds, currentObservation );


////    Eigen::VectorXd finalRelativePositionIoEarth = variationalEquationsSimulator.getDynamicsSimulator( )
////            ->getDependentVariableHistory( ).rbegin( )->second.segment( 0, 3 );
////    Eigen::VectorXd centralInstantRelativePositionEuropaEarth = variationalEquationsSimulator.getDynamicsSimulator( )
////            ->getDependentVariableHistory( ).rbegin( )->second.segment( 3, 3 );


//    // Right ascensions and declinations Io
//    Eigen::Matrix< double, 3, 1 > prepreFinalSphericalStateIo =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//               prepreFinalRelativeStateIoEarth.segment( 0, 3 ) ).template cast< double >( );
//    double prepreFinalRightAscensionIo = prepreFinalSphericalStateIo.z( );
//    double prepreFinalDeclinationIo = mathematical_constants::PI / 2.0 - prepreFinalSphericalStateIo.y( );

//    Eigen::Matrix< double, 3, 1 > preFinalSphericalStateIo =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                preFinalRelativeStateIoEarth.segment( 0, 3 ) ).template cast< double >( );
//    double preFinalRightAscensionIo = preFinalSphericalStateIo.z( );
//    double preFinalDeclinationIo = mathematical_constants::PI / 2.0 - preFinalSphericalStateIo.y( );

//    Eigen::Matrix< double, 3, 1 > finalSphericalStateIo =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                finalRelativeStateIoEarth.segment( 0, 3 ) ).template cast< double >( );
//    double finalRightAscensionIo = finalSphericalStateIo.z( );
//    double finalDeclinationIo = mathematical_constants::PI / 2.0 - finalSphericalStateIo.y( );


//    // Right ascensions and declinations Europa
//    Eigen::Matrix< double, 3, 1 > prepreFinalSphericalStateEuropa =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//               prepreFinalRelativeStateEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
//    double prepreFinalRightAscensionEuropa = prepreFinalSphericalStateEuropa.z( );
//    double prepreFinalDeclinationEuropa = mathematical_constants::PI / 2.0 - prepreFinalSphericalStateEuropa.y( );

//    Eigen::Matrix< double, 3, 1 > preFinalSphericalStateEuropa =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                preFinalRelativeStateEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
//    double preFinalRightAscensionEuropa = preFinalSphericalStateEuropa.z( );
//    double preFinalDeclinationEuropa = mathematical_constants::PI / 2.0 - preFinalSphericalStateEuropa.y( );

//    Eigen::Matrix< double, 3, 1 > finalSphericalStateEuropa =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                finalRelativeStateEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
//    double finalRightAscensionEuropa = finalSphericalStateEuropa.z( );
//    double finalDeclinationEuropa = mathematical_constants::PI / 2.0 - finalSphericalStateEuropa.y( );



//    // Variations in right ascension and declination for Io
//    double variationRightAscension = finalRightAscensionIo - preFinalRightAscensionIo;
//    double variationDeclination = finalDeclinationIo - preFinalDeclinationIo;

//    double preVariationRightAscension = preFinalRightAscensionIo - prepreFinalRightAscensionIo;
//    double preVariationDeclination = preFinalDeclinationIo - prepreFinalDeclinationIo;

////    double variationX = computeX(  );

//    Eigen::Vector3d finalRelativePositionIoEarth = finalRelativeStateIoEarth.segment( 0, 3 );
//    Eigen::Vector3d finalRelativeVelocityIoEarth = finalRelativeStateIoEarth.segment( 3, 3 );
//    double partialRightAscension = computePartialRightAscensionWrtTime(
//                finalRelativePositionIoEarth, finalRelativeVelocityIoEarth );

//    double partialDeclination = computePartialDeclinationWrtTime(
//                finalRelativePositionIoEarth, finalRelativeVelocityIoEarth );

//    Eigen::Vector3d preFinalRelativePositionIoEarth = preFinalRelativeStateIoEarth.segment( 0, 3 );
//    Eigen::Vector3d preFinalRelativeVelocityIoEarth = preFinalRelativeStateIoEarth.segment( 3, 3 );
//    double prePartialRightAscension = computePartialRightAscensionWrtTime(
//                preFinalRelativePositionIoEarth, preFinalRelativeVelocityIoEarth );

//    double prePartialDeclination = computePartialDeclinationWrtTime(
//                preFinalRelativePositionIoEarth, preFinalRelativeVelocityIoEarth );

//    double approximatePartialRightAscension = variationRightAscension / timeVariation;
//    double approximatePartialDeclination = variationDeclination / timeVariation;

//    std::cout << "variation right ascension: " << variationRightAscension << "\n\n";
//    std::cout << "approximate variation right ascension: " << timeVariation * prePartialRightAscension << "\n\n";

//    std::cout << "variation declination: " << variationDeclination << "\n\n";
//    std::cout << "approximate variation declination: " << timeVariation * prePartialDeclination << "\n\n";

//    std::cout << "partial right ascension: " << prePartialRightAscension << "\n\n";
//    std::cout << "approximate partial right ascension: " << approximatePartialRightAscension << "\n\n";

//    std::cout << "partial declination: " << prePartialDeclination << "\n\n";
//    std::cout << "approximate partial declination: " << approximatePartialDeclination << "\n\n";


//    Eigen::Vector3d prepreFinalRelativePositionIoEarth = prepreFinalRelativeStateIoEarth.segment( 0, 3 );
//    Eigen::Vector3d prepreFinalRelativeVelocityIoEarth = prepreFinalRelativeStateIoEarth.segment( 3, 3 );
//    double preprePartialRightAscension = computePartialRightAscensionWrtTime(
//                prepreFinalRelativePositionIoEarth, prepreFinalRelativeVelocityIoEarth );

//    double preprePartialDeclination = computePartialDeclinationWrtTime(
//                prepreFinalRelativePositionIoEarth, prepreFinalRelativeVelocityIoEarth );

//    double preApproximatePartialRightAscension = preVariationRightAscension / preTimeVariation;
//    double preApproximatePartialDeclination = preVariationDeclination / preTimeVariation;

////    std::cout << "variation right ascension: " << preVariationRightAscension << "\n\n";
////    std::cout << "approximate variation right ascension: " << preTimeVariation * preprePartialRightAscension << "\n\n";

////    std::cout << "variation declination: " << preVariationDeclination << "\n\n";
////    std::cout << "approximate variation declination: " << preTimeVariation * preprePartialDeclination << "\n\n";

////    std::cout << "partial right ascension: " << preprePartialRightAscension << "\n\n";
////    std::cout << "approximate partial right ascension: " << preApproximatePartialRightAscension << "\n\n";

////    std::cout << "partial declination: " << preprePartialDeclination << "\n\n";
////    std::cout << "approximate partial declination: " << preApproximatePartialDeclination << "\n\n";

//    Eigen::Vector3d preFinalRelativePositionEuropaEarth = preFinalRelativeStateEuropaEarth.segment( 0, 3 );
//    Eigen::Vector3d preFinalRelativeVelocityEuropaEarth = preFinalRelativeStateEuropaEarth.segment( 3, 3 );

//    double variationX = ( computeX( finalRightAscensionIo, finalDeclinationIo, finalRightAscensionEuropa, finalDeclinationEuropa )
//            - computeX( preFinalRightAscensionIo, preFinalDeclinationIo, preFinalRightAscensionEuropa, preFinalDeclinationEuropa ) ) * 3600.0 * 180.0 / mathematical_constants::PI;
//    double variationY = ( computeY( finalDeclinationIo, finalDeclinationEuropa ) - computeY( preFinalDeclinationIo, preFinalDeclinationEuropa ) ) * 3600.0 * 180.0 / mathematical_constants::PI;

//    double approximatePartialX = variationX / timeVariation;
//    double approximatePartialY = variationY / timeVariation;

//    double partialX = computePartialXwrtTime( preFinalRelativePositionIoEarth, preFinalRelativePositionEuropaEarth, preFinalRelativeVelocityIoEarth, preFinalRelativeVelocityEuropaEarth,
//                                              preFinalRightAscensionIo, preFinalDeclinationIo, preFinalRightAscensionEuropa, preFinalDeclinationEuropa ) * 3600.0 * 180.0 / mathematical_constants::PI;
//    double partialY = computePartialYwrtTime( preFinalRelativePositionIoEarth, preFinalRelativePositionEuropaEarth, preFinalRelativeVelocityIoEarth, preFinalRelativeVelocityEuropaEarth ) * 3600.0 * 180.0 / mathematical_constants::PI;

//    std::cout << "value X: " << computeX( preFinalRightAscensionIo, preFinalDeclinationIo, preFinalRightAscensionEuropa, preFinalDeclinationEuropa ) << "\n\n";
//    std::cout << "value Y: " << computeY( preFinalDeclinationIo, preFinalDeclinationEuropa ) << "\n\n";

//    std::cout << "partial X: " << partialX << "\n\n";
//    std::cout << "approximate partial X: " << approximatePartialX << "\n\n";
//    std::cout << "partial Y: " << partialY << "\n\n";
//    std::cout << "approximate partial Y: " << approximatePartialY << "\n\n";



//    /// TEST SECOND DERIVATIVES

//    double variationPartialRightAscension = //prePartialRightAscension - preprePartialRightAscension; //approximatePartialRightAscension - preApproximatePartialRightAscension;
//            partialRightAscension - preprePartialRightAscension;

//    double approximateSecondPartialRightAscension = variationPartialRightAscension / ( 2.0 * timeVariation );

//    double variationPartialDeclination = //prePartialDeclination - preprePartialDeclination;// approximatePartialDeclination - preApproximatePartialDeclination;
//            partialDeclination - preprePartialDeclination;
//    double approximateSecondPartialDeclination = variationPartialDeclination / ( 2.0 * timeVariation);

//    double secondPartialRightAscension = computeSecondPartialRightAscensionWrtTime( preFinalRelativePositionIoEarth,
//                                                                                    preFinalRelativeVelocityIoEarth,
//                                                                                    preFinalRelativeAccelerationIoEarth );

//    double secondPartialDeclination = computeSecondPartialDeclinationWrtTime( preFinalRelativePositionIoEarth,
//                                                                              preFinalRelativeVelocityIoEarth,
//                                                                              preFinalRelativeAccelerationIoEarth );

//    std::cout << "second partial right ascension: " << secondPartialRightAscension << "\n\n";
//    std::cout << "approximate partial right ascension: " << approximateSecondPartialRightAscension << "\n\n";

//    std::cout << "second partial declination: " << secondPartialDeclination << "\n\n";
//    std::cout << "approximate partial declination: " << approximateSecondPartialDeclination << "\n\n";





    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        TEST PARTIAL RIGHT ASCENSION/DECLINATION WRT POSITION      /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    double estimatedCentralInstant = 1341680.0; //1342800.0; //1344949.7809586974327; //1342800.0; // -72932.32641653056 -48605.27593299343 + 6494.533911467588 -71.65702315459203 -0.01173150924660149;


//    // Create light-time correction settings
//    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
//    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
//    lightTimeCorrectionSettings.push_back( std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
//                                                lightTimePerturbingBodies ) );

//    // Create observation settings
//    std::shared_ptr< observation_models::ObservationSettings > observableSettings = std::make_shared< observation_models::ObservationSettings >
//            ( observation_models::mutual_approximation, lightTimeCorrectionSettings, std::make_shared< observation_models::ConstantObservationBiasSettings >(
//                  ( Eigen::Vector1d( ) << 0.0 ).finished( ), true ) );

//    // Create observation model.
//    std::shared_ptr< observation_models::ObservationModel< 1, double, double > > observationModel =
//           observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
//                linkEnds, observableSettings, bodyMap );

//    std::vector< double > linkEndTimes;
//    linkEndStates.clear( );
//    Eigen::Vector1d observationFromReceptionTime = observationModel->computeObservationsWithLinkEndData( estimatedCentralInstant, observation_models::receiver, linkEndTimes, linkEndStates );
//    std::cout << "central instant from observation model: " << observationFromReceptionTime << "\n\n";
//    std::cout << "estimated central instant after iterations: " << 1342800.0 -72932.32641653056 -48605.27593299343 + 6494.533911467588 -71.65702315459203 -0.01173150924660149 << "\n\n";

//    std::cout << "LINK ENDS TIMES: " << linkEndTimes[ 0 ] << " & " << linkEndTimes[ 1 ] << " & " << linkEndTimes[ 2 ] << "\n\n";

    std::cout << "state transition matrix at estimated central instant: " << stateTransitionInterpolator->interpolate( estimatedCentralInstant ) << "\n\n";

    Eigen::VectorXd centralInstantState = stateHistoryInterpolator->interpolate( estimatedCentralInstant );
    Eigen::VectorXd centralInstantDependentVariables = dependentVariablesInterpolator->interpolate( estimatedCentralInstant );
    Eigen::Vector6d centralInstantEarthState = centralInstantDependentVariables.segment( 186, 6 ); // currentPropagatedState.segment( 0, 6 ); //  systemInitialState.segment( 0, 6 );
    Eigen::Vector6d centralInstantIoState = centralInstantDependentVariables.segment( 198, 6 ); //currentPropagatedState.segment( 12, 6 ); //systemInitialState.segment( 12, 6 );
    Eigen::Vector6d centralInstantEuropaState = centralInstantDependentVariables.segment( 204, 6 ); //currentPropagatedState.segment( 18, 6 ); //systemInitialState.segment( 18, 6 );

//    std::cout << "TEST SYSTEM INITIAL STATE: " << systemInitialState.transpose( ) -  stateHistoryInterpolator->interpolate( simulationStartEpoch ).transpose( ) << "\n\n";
//    std::cout << "TEST: " << stateHistoryInterpolator->interpolate( simulationStartEpoch ).transpose( ) << "\n\n";

//    std::pair< double, double > currentRightAscensionDeclinationIo = computeRightAscensionDeclination( centralInstantDependentVariables.segment( 9, 3 ) );
//    std::pair< double, double > currentRightAscensionDeclinationEuropa = computeRightAscensionDeclination( centralInstantDependentVariables.segment( 15, 3 ) );
//    double currentAverageDeclination = ( currentRightAscensionDeclinationIo.second + currentRightAscensionDeclinationEuropa.second ) / 2.0;
//    std::cout << "apparent distance Io-Europa (arcseconds): " << std::sqrt(
//                     ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                     * ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                     + ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second )
//                     * ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second ) ) * 3600.0 * 180.0 / mathematical_constants::PI << " at time = " << itr->first << " seconds." "\n\n";
//    currentRightAscensionDeclinationIo = computeRightAscensionDeclination( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) );
//    currentRightAscensionDeclinationEuropa = computeRightAscensionDeclination( ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
//    currentAverageDeclination = ( currentRightAscensionDeclinationIo.second + currentRightAscensionDeclinationEuropa.second ) / 2.0;
//    std::cout << "VALIDATION apparent distance Io-Europa (arcseconds): " << std::sqrt(
//                     ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                     * ( currentRightAscensionDeclinationEuropa.first - currentRightAscensionDeclinationIo.first ) * std::cos( currentAverageDeclination )
//                     + ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second )
//                     * ( currentRightAscensionDeclinationEuropa.second - currentRightAscensionDeclinationIo.second ) ) * 3600.0 * 180.0 / mathematical_constants::PI << " at time = " << itr->first << " seconds." "\n\n";

    // Update mutual approximation scaling.
//    linkEndStates.clear( );

    observation_models::LinkEnds linkEnds;
    linkEnds[ observation_models::receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ observation_models::transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ observation_models::transmitter2 ] = std::make_pair( "Europa" , ""  );

    std::vector< Eigen::Vector6d > linkEndStates;
    linkEndStates.push_back( centralInstantIoState );
    linkEndStates.push_back( centralInstantEuropaState );
    linkEndStates.push_back( centralInstantEarthState );

//    times.clear( );
    std::vector< double > linkEndsTimes;
    linkEndsTimes.push_back( estimatedCentralInstant /*simulationStartEpoch*/ );
    linkEndsTimes.push_back( estimatedCentralInstant /*simulationStartEpoch*/ );
    linkEndsTimes.push_back( estimatedCentralInstant /*simulationStartEpoch*/ );

    Eigen::VectorXd currentObservation;
    mutualApproximationScaling->update( linkEndStates, linkEndsTimes, observation_models::receiver, linkEnds, currentObservation );



    // Create light-time correction settings
    std::vector< std::string > perturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                perturbingBodies ) );

    // Create observation settings
    std::shared_ptr< observation_models::ObservationSettings > observableSettings = std::make_shared< observation_models::ObservationSettings >
            ( observation_models::apparent_distance, lightTimeCorrectionSettings, std::make_shared< observation_models::ConstantObservationBiasSettings >(
                  ( Eigen::Vector1d( ) << 0.0 ).finished( ), true ) );

    // Create apparent distance observation model.
    std::shared_ptr< observation_models::ObservationModel< 1, double, double > > apparentDistanceModel =
           observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, observableSettings, bodyMap );

    std::map< double, Eigen::Vector1d > apparentDistanceHistory;
    for ( std::map< double, Eigen::VectorXd >::iterator itr = stateHistory.begin( ) ; itr!= stateHistory.end( ) ; itr++ )
    {
        Eigen::Vector3d positionIoEarth = dependentVariablesHistory[ itr->first ].segment( 9, 3 );
        Eigen::Vector3d positionEuropaEarth = dependentVariablesHistory[ itr->first ].segment( 15, 3 );
        std::pair< double, double > angularPositionIoWrtEarth = computeRightAscensionAndDeclination( positionIoEarth );
        std::pair< double, double > angularPositionEuropaWrtEarth = computeRightAscensionAndDeclination( positionEuropaEarth );
        apparentDistanceHistory[ itr->first ] = (
                    Eigen::Vector1d( ) << 3600.0 * 180.0 / mathematical_constants::PI * /*apparentDistanceModel->computeObservations( itr->first, observation_models::receiver )*/
                    std::sqrt( ( angularPositionEuropaWrtEarth.first - angularPositionIoWrtEarth.first )
                               * std::cos( ( angularPositionEuropaWrtEarth.second + angularPositionIoWrtEarth.second ) / 2.0 )
                               * ( angularPositionEuropaWrtEarth.first - angularPositionIoWrtEarth.first )
                               * std::cos( ( angularPositionEuropaWrtEarth.second + angularPositionIoWrtEarth.second ) / 2.0 )
                               + ( angularPositionEuropaWrtEarth.second - angularPositionIoWrtEarth.second )
                               * ( angularPositionEuropaWrtEarth.second - angularPositionIoWrtEarth.second ) ) ).finished( );
    }


    input_output::writeDataMapToTextFile( apparentDistanceHistory,
                                          "apparentDistanceHistory.dat",
                                          "C:/Users/chamb/Documents/PhD/testMutualApproximationPartials/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );



//    Eigen::Matrix< double, 3, 1 > initialSphericalStateIo =
//            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) ).template cast< double >( );

    std::pair< double, double > rightAscensionAndDeclinationIo = computeRightAscensionDeclination( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) );
    std::pair< double, double > rightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );

//    std::cout << "initial right ascension Io: " << rightAscensionAndDeclinationIo.first << "\n\n";
//    std::cout << "initial right ascension Io - validation: " << initialSphericalStateIo.z( ) << "\n\n";

//    std::cout << "initial declination Io: " << rightAscensionAndDeclinationIo.second << "\n\n";
//    std::cout << "initial declination Io - validation: " << mathematical_constants::PI / 2.0 - initialSphericalStateIo.y( ) << "\n\n";

    Eigen::Vector2d instrumentalFrameRelativePosition = mutualApproximationScaling->getInstrumentalFrameRelativePosition( );
    double centralInstantValueX = instrumentalFrameRelativePosition[ 0 ]; //computeX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                     ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
    double centralInstantValueY = instrumentalFrameRelativePosition[ 1 ]; //computeY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                     ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );

    Eigen::Vector3d partialRightAscensionWrtPositionIo = computePartialOfRightAscensionWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), true );
    Eigen::Vector3d partialDeclinationWrtPositionIo = computePartialOfDeclinationWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), true );

    Eigen::Vector3d positionPartialOfPartialRightAscensionWrtTime =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    Eigen::Vector3d positionPartialOfPartialDeclinationWrtTime =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                              ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );


//    Eigen::Vector3d partialOfXwrtPosition = computePositionPartialOfX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                   ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ), true );
    Eigen::Matrix< double, 2, 3 > partialOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition
            = mutualApproximationScaling->getPartialsOfRelativePositionInInstrumentalFrameWrtFirstTransmitterPosition( );
    Eigen::Vector3d partialOfXwrtPosition = partialOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition.block( 0, 0, 1, 3 ).transpose( );
//    Eigen::Vector3d partialOfYwrtPosition = computePositionPartialOfY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                   ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ), true );
    Eigen::Vector3d partialOfYwrtPosition = partialOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition.block( 1, 0, 1, 3 ).transpose( );


    Eigen::Vector2d instrumentalFrameRelativeVelocity = mutualApproximationScaling->getInstrumentalFrameRelativeVelocity( );
//    double firstTimePartialOfX = computePartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                 rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                                                 rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    double firstTimePartialOfX = instrumentalFrameRelativeVelocity[ 0 ];
//    double firstTimePartialOfY = computePartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );
    double firstTimePartialOfY = instrumentalFrameRelativeVelocity[ 1 ];


//    std::cout << "value X: " << computeX( rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                          rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second ) << "\n\n";
//    std::cout << "value Y: " << computeY( rightAscensionAndDeclinationIo.second, rightAscensionAndDeclinationEuropa.second ) << "\n\n";
//    std::cout << "value velocity Vx in instrumental frame: " << firstTimePartialOfX << "\n\n";
//    std::cout << "value velocity Vy in instrumental frame: " << firstTimePartialOfY << "\n\n";
//    std::cout << "value acceleration Ax in instrumental frame: " <<
//                 computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                               ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                               centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 0, 3 ) - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 6, 3 ),
//                                               centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 3, 3 ) - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 6, 3 ),
//                                               rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                               rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second ) << "\n\n";
//    std::cout << "value acceleration Ay in instrumental frame: " <<
//                 computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                               ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                               centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 0, 3 )
//                                               - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 6, 3 ),
//                                               centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 3, 3 )
//                                               - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 6, 3 ) ) << "\n\n";


    Eigen::Matrix< double, 2, 3 > partialOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition
            = mutualApproximationScaling->getPartialsOfRelativeVelocityInInstrumentalFrameWrtFirstTransmitterPosition( );
    Eigen::Vector3d positionPartialOfFirstTimeDerivativeOfX = partialOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition.block( 0, 0, 1, 3 ).transpose( );
    Eigen::Vector3d positionPartialOfFirstTimeDerivativeOfY = partialOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition.block( 1, 0, 1, 3 ).transpose( );
//    Eigen::Vector3d positionPartialOfFirstTimeDerivativeOfX = computepositionPartialOfFirstTimeDerivativeOfX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                                                          ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                                                          ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                                                          ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ), true );
//    Eigen::Vector3d positionPartialOfFirstTimeDerivativeOfY = computepositionPartialOfFirstTimeDerivativeOfY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                                                          ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                                                          ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                                                          ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ), true );



    Eigen::Vector3d variationCentralInstantPositionIoEarth = 0.0001 * centralInstantState.segment( 12, 3 ); // 1.0e-7 * ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 );
    Eigen::Vector3d variationCentralInstantPositionEuropaEarth = 0.0001 * centralInstantState.segment( 18, 3 ); //1.0e-7 * ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 );




    /// x - axis
    Eigen::Vector3d modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +
            ( Eigen::Vector3d( ) << variationCentralInstantPositionIoEarth[ 0 ], 0.0, 0.0 ).finished( );
    std::pair< double, double > modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );


    std::vector< Eigen::Vector6d > modifiedLinkEndStates;
    Eigen::Vector6d modifiedCentralInstantIo = centralInstantIoState;
    modifiedCentralInstantIo[ 0 ] += variationCentralInstantPositionIoEarth[ 0 ];
    modifiedLinkEndStates.push_back( modifiedCentralInstantIo );
    modifiedLinkEndStates.push_back( centralInstantEuropaState );
    modifiedLinkEndStates.push_back( centralInstantEarthState );

    Eigen::VectorXd modifiedObservation;
    modifiedMutualApproximationScaling->update( modifiedLinkEndStates, linkEndsTimes, observation_models::receiver, linkEnds, modifiedObservation );


    double variationRightAscensionWrtPosition = ( modifiedRightAscensionAndDeclinationIo.first - rightAscensionAndDeclinationIo.first )
            / variationCentralInstantPositionIoEarth[ 0 ];
    double variationDeclinationWrtPosition = ( modifiedRightAscensionAndDeclinationIo.second - rightAscensionAndDeclinationIo.second )
            / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical partial right ascension wrt position - x: " << variationRightAscensionWrtPosition << "\n\n";
    std::cout << "validation: " << partialRightAscensionWrtPositionIo[ 0 ] << "\n\n";
    std::cout << "numerical partial declination wrt position - x: " << variationDeclinationWrtPosition << "\n\n";
    std::cout << "validation: " << partialDeclinationWrtPositionIo[ 0 ] << "\n\n";


    double partialRightAscensionWrtTime = computePartialOfRightAscensionWrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                               ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    double modifiedPartialRightAscensionWrtTime = computePartialOfRightAscensionWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                       ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    double variationPartialRightAscensionWrtTime = ( modifiedPartialRightAscensionWrtTime - partialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 0 ];

    double partialDeclinationWrtTime = computePartialOfDeclinationWrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                         ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    double modifiedPartialDeclinationWrtTime = computePartialOfDeclinationWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    double variationPartialDeclinationWrtTime = ( modifiedPartialDeclinationWrtTime - partialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of first partial right ascension wrt time - x: " << variationPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialRightAscensionWrtTime[ 0 ] << "\n\n";
    std::cout << "numerical position partial of first partial declination wrt time - x: " << variationPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialDeclinationWrtTime[ 0 ] << "\n\n";


    double modifiedCentralInstantValueX = ( rightAscensionAndDeclinationEuropa.first - modifiedRightAscensionAndDeclinationIo.first )
            * std::cos( ( rightAscensionAndDeclinationEuropa.second + modifiedRightAscensionAndDeclinationIo.second ) / 2.0 );
//    instrumentalFrameRelativePosition_[ 1 ] = declinationSecondTransmitter_ - declinationFirstTransmitter_;
    double modifiedcentralInstantValueY = rightAscensionAndDeclinationEuropa.second - modifiedRightAscensionAndDeclinationIo.second;

//    double modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
    double variationValueXwrtIoPosition = ( modifiedCentralInstantValueX - centralInstantValueX ) / variationCentralInstantPositionIoEarth[ 0 ];
//    double modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
    double variationValueYwrtIoPosition = ( modifiedcentralInstantValueY - centralInstantValueY ) / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of X - x: " << variationValueXwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfXwrtPosition[ 0 ] << "\n\n";
    std::cout << "numerical position partial of Y - x: " << variationValueYwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfYwrtPosition[ 0 ] << "\n\n";

    Eigen::Vector2d modifiedInstrumentalFrameRelativeVelocity = modifiedMutualApproximationScaling->getInstrumentalFrameRelativeVelocity( );
    double modifiedFirstTimeDerivativeOfX = modifiedInstrumentalFrameRelativeVelocity[ 0 ];
//            computePartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                    ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                    modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                    rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    double variationFirstTimePartiaOfX = ( modifiedFirstTimeDerivativeOfX - firstTimePartialOfX ) / variationCentralInstantPositionIoEarth[ 0 ];

    double modifiedFirstTimeDerivativeOfY = modifiedInstrumentalFrameRelativeVelocity[ 1 ]; //computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );
    double variationFirstTimePartiaOfY = ( modifiedFirstTimeDerivativeOfY - firstTimePartialOfY ) / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of first partial X wrt time - x: " << variationFirstTimePartiaOfX << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfX[ 0 ] << "\n\n";
    std::cout << "numerical position partial of first partial Y wrt time - x: " << variationFirstTimePartiaOfY << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfY[ 0 ] << "\n\n";



    /// y - axis
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +
            ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionIoEarth[ 1 ], 0.0 ).finished( );
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );

    modifiedLinkEndStates.clear( );
    modifiedCentralInstantIo = centralInstantIoState;
    modifiedCentralInstantIo[ 1 ] += variationCentralInstantPositionIoEarth[ 1 ];
    modifiedLinkEndStates.push_back( modifiedCentralInstantIo );
    modifiedLinkEndStates.push_back( centralInstantEuropaState );
    modifiedLinkEndStates.push_back( centralInstantEarthState );
    modifiedMutualApproximationScaling->update( modifiedLinkEndStates, linkEndsTimes, observation_models::receiver, linkEnds, modifiedObservation );

    variationRightAscensionWrtPosition = ( modifiedRightAscensionAndDeclinationIo.first - rightAscensionAndDeclinationIo.first )
            / variationCentralInstantPositionIoEarth[ 1 ];
    variationDeclinationWrtPosition = ( modifiedRightAscensionAndDeclinationIo.second - rightAscensionAndDeclinationIo.second )
            / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical partial right ascension wrt position - y: " << variationRightAscensionWrtPosition << "\n\n";
    std::cout << "validation: " << partialRightAscensionWrtPositionIo[ 1 ] << "\n\n";
    std::cout << "numerical partial declination wrt position - y: " << variationDeclinationWrtPosition << "\n\n";
    std::cout << "validation: " << partialDeclinationWrtPositionIo[ 1 ] << "\n\n";

    modifiedPartialRightAscensionWrtTime = computePartialOfRightAscensionWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                       ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    variationPartialRightAscensionWrtTime = ( modifiedPartialRightAscensionWrtTime - partialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 1 ];

    modifiedPartialDeclinationWrtTime = computePartialOfDeclinationWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    variationPartialDeclinationWrtTime = ( modifiedPartialDeclinationWrtTime - partialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of first partial right ascension wrt time - y: " << variationPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialRightAscensionWrtTime[ 1 ] << "\n\n";
    std::cout << "numerical position partial of first partial declination wrt time - y: " << variationPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialDeclinationWrtTime[ 1 ] << "\n\n";


//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
    modifiedCentralInstantValueX = ( rightAscensionAndDeclinationEuropa.first - modifiedRightAscensionAndDeclinationIo.first )
            * std::cos( ( rightAscensionAndDeclinationEuropa.second + modifiedRightAscensionAndDeclinationIo.second ) / 2.0 );
    variationValueXwrtIoPosition = ( modifiedCentralInstantValueX - centralInstantValueX ) / variationCentralInstantPositionIoEarth[ 1 ];
    modifiedcentralInstantValueY = rightAscensionAndDeclinationEuropa.second - modifiedRightAscensionAndDeclinationIo.second;
    variationValueYwrtIoPosition = ( modifiedcentralInstantValueY - centralInstantValueY ) / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of X - y: " << variationValueXwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfXwrtPosition[ 1 ] << "\n\n";
    std::cout << "numerical position partial of Y - y: " << variationValueYwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfYwrtPosition[ 1 ] << "\n\n";


    modifiedInstrumentalFrameRelativeVelocity = modifiedMutualApproximationScaling->getInstrumentalFrameRelativeVelocity( );
    modifiedFirstTimeDerivativeOfX = modifiedInstrumentalFrameRelativeVelocity[ 0 ];
//    modifiedFirstTimeDerivativeOfX =
//            computePartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                    ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                    modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                    rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    variationFirstTimePartiaOfX = ( modifiedFirstTimeDerivativeOfX - firstTimePartialOfX ) / variationCentralInstantPositionIoEarth[ 1 ];

    modifiedFirstTimeDerivativeOfY = modifiedInstrumentalFrameRelativeVelocity[ 1 ];
//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );
    variationFirstTimePartiaOfY = ( modifiedFirstTimeDerivativeOfY - firstTimePartialOfY ) / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of first partial X wrt time - y: " << variationFirstTimePartiaOfX << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfX[ 1 ] << "\n\n";
    std::cout << "numerical position partial of first partial Y wrt time - y: " << variationFirstTimePartiaOfY << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfY[ 1 ] << "\n\n";




    /// z - axis
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +
            ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionIoEarth[ 2 ] ).finished( );
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );

    modifiedLinkEndStates.clear( );
    modifiedCentralInstantIo = centralInstantIoState;
    modifiedCentralInstantIo[ 2 ] += variationCentralInstantPositionIoEarth[ 2 ];
    modifiedLinkEndStates.push_back( modifiedCentralInstantIo );
    modifiedLinkEndStates.push_back( centralInstantEuropaState );
    modifiedLinkEndStates.push_back( centralInstantEarthState );
    modifiedMutualApproximationScaling->update( modifiedLinkEndStates, linkEndsTimes, observation_models::receiver, linkEnds, modifiedObservation );

    variationRightAscensionWrtPosition = ( modifiedRightAscensionAndDeclinationIo.first - rightAscensionAndDeclinationIo.first )
            / variationCentralInstantPositionIoEarth[ 2 ];
    variationDeclinationWrtPosition = ( modifiedRightAscensionAndDeclinationIo.second - rightAscensionAndDeclinationIo.second )
            / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical partial right ascension wrt position - z: " << variationRightAscensionWrtPosition << "\n\n";
    std::cout << "validation: " << partialRightAscensionWrtPositionIo[ 2 ] << "\n\n";
    std::cout << "numerical partial declination wrt position - z: " << variationDeclinationWrtPosition << "\n\n";
    std::cout << "validation: " << partialDeclinationWrtPositionIo[ 2 ] << "\n\n";


    modifiedPartialRightAscensionWrtTime = computePartialOfRightAscensionWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                       ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    variationPartialRightAscensionWrtTime = ( modifiedPartialRightAscensionWrtTime - partialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 2 ];

    modifiedPartialDeclinationWrtTime = computePartialOfDeclinationWrtTime( modifiedCentralInstantPositionIoEarth,
                                                                                 ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ) );
    variationPartialDeclinationWrtTime = ( modifiedPartialDeclinationWrtTime - partialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of first partial right ascension wrt time - z: " << variationPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialRightAscensionWrtTime[ 2 ] << "\n\n";
    std::cout << "numerical position partial of first partial declination wrt time - z: " << variationPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << positionPartialOfPartialDeclinationWrtTime[ 2 ] << "\n\n";


    modifiedCentralInstantValueX = ( rightAscensionAndDeclinationEuropa.first - modifiedRightAscensionAndDeclinationIo.first )
            * std::cos( ( rightAscensionAndDeclinationEuropa.second + modifiedRightAscensionAndDeclinationIo.second ) / 2.0 );
    variationValueXwrtIoPosition = ( modifiedCentralInstantValueX - centralInstantValueX ) / variationCentralInstantPositionIoEarth[ 2 ];
    modifiedcentralInstantValueY = rightAscensionAndDeclinationEuropa.second - modifiedRightAscensionAndDeclinationIo.second;
    variationValueYwrtIoPosition = ( modifiedcentralInstantValueY - centralInstantValueY ) / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of X - z: " << variationValueXwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfXwrtPosition[ 2 ] << "\n\n";
    std::cout << "numerical position partial of Y - z: " << variationValueYwrtIoPosition << "\n\n";
    std::cout << "validation: " << partialOfYwrtPosition[ 2 ] << "\n\n";


    modifiedInstrumentalFrameRelativeVelocity = modifiedMutualApproximationScaling->getInstrumentalFrameRelativeVelocity( );
    modifiedFirstTimeDerivativeOfX = modifiedInstrumentalFrameRelativeVelocity[ 0 ];
//    modifiedFirstTimeDerivativeOfX =
//            computePartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                    ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                    modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                    rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    variationFirstTimePartiaOfX = ( modifiedFirstTimeDerivativeOfX - firstTimePartialOfX ) / variationCentralInstantPositionIoEarth[ 2 ];

    modifiedFirstTimeDerivativeOfY = modifiedInstrumentalFrameRelativeVelocity[ 1 ];
//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );
    variationFirstTimePartiaOfY = ( modifiedFirstTimeDerivativeOfY - firstTimePartialOfY ) / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of first partial X wrt time - z: " << variationFirstTimePartiaOfX << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfX[ 2 ] << "\n\n";
    std::cout << "numerical position partial of first partial Y wrt time - z: " << variationFirstTimePartiaOfY << "\n\n";
    std::cout << "validation: " << positionPartialOfFirstTimeDerivativeOfY[ 2 ] << "\n\n";




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////        TEST POSITION PARTIAL OF SECOND PARTIAL RIGHT ASCENSION/DECLINATION WRT TIME      /////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::string > sunCentralBodies;
    sunCentralBodies.push_back( "Sun" );
    sunCentralBodies.push_back( "Sun" );
    sunCentralBodies.push_back( "Sun" );
    sunCentralBodies.push_back( "Sun" );

    Eigen::VectorXd sunCenteredState = centralInstantDependentVariables.segment( 186, 24 ); // propagators::getInitialStatesOfBodies(
//                bodiesToPropagate, sunCentralBodies, bodyMap, simulationStartEpoch );

//    std::cout << "difference initial state: " << ( sunCenteredState - systemInitialState ).transpose( ) << "\n\n";

    Eigen::Vector3d relativePositionJupiterSun = - sunCenteredState.segment( 6, 3 );
    Eigen::Vector3d relativePositionJupiterIo = sunCenteredState.segment( 12, 3 ) - sunCenteredState.segment( 6, 3 );
    Eigen::Vector3d relativePositionJupiterEuropa = sunCenteredState.segment( 18, 3 ) - sunCenteredState.segment( 6, 3 );

    Eigen::Vector3d centralInstantAccelerationJupiter = centralInstantDependentVariables.segment( 183, 3 ); //Eigen::Vector3d::Zero( );
//    centralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationIo = dependentVariables.begin( )->second.segment( 0, 3 ); //.transpose( );

//    std::cout << "initial acceleration Jupiter: " << centralInstantAccelerationJupiter.transpose( ) << "\n\n";
//    std::cout << "validation: " << centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 147, 3 ).transpose( ) << "\n\n";

    Eigen::Vector3d relativePositionIoSun = - sunCenteredState.segment( 12, 3 );
    Eigen::Vector3d relativePositionIoEarth = sunCenteredState.segment( 0, 3 ) - sunCenteredState.segment( 12, 3 );
    Eigen::Vector3d relativePositionIoJupiter = sunCenteredState.segment( 6, 3 ) - sunCenteredState.segment( 12, 3 );
    Eigen::Vector3d relativePositionIoEuropa = sunCenteredState.segment( 18, 3 ) - sunCenteredState.segment( 12, 3 );

    Eigen::Vector3d centralInstantAccelerationIo = centralInstantDependentVariables.segment( 0, 3 ); //Eigen::Vector3d::Zero( );
//    centralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationIo = ( centralInstantAccelerationIo - centralInstantAccelerationJupiter ); //dependentVariables.begin( )->second.segment( 0, 3 ); //.transpose( );

//    std::cout << "initial acceleration Io: " << centralInstantAccelerationIo.transpose( ) << "\n\n";
//    std::cout << "validation: " << centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 0, 3 ).transpose( ) << "\n\n";
//    std::cout << "test initial acceleration Io: " << ( centralInstantAccelerationIo - centralInstantAccelerationJupiter ).transpose( ) << "\n\n";

    Eigen::Vector3d relativePositionEuropaSun = - sunCenteredState.segment( 18, 3 );
    Eigen::Vector3d relativePositionEuropaEarth = sunCenteredState.segment( 0, 3 ) - sunCenteredState.segment( 18, 3 );
    Eigen::Vector3d relativePositionEuropaJupiter = sunCenteredState.segment( 6, 3 ) - sunCenteredState.segment( 18, 3 );
    Eigen::Vector3d relativePositionEuropaIo = sunCenteredState.segment( 12, 3 ) - sunCenteredState.segment( 18, 3 );

    Eigen::Vector3d centralInstantAccelerationEuropa = centralInstantDependentVariables.segment( 3, 3 ); //Eigen::Vector3d::Zero( );
//    centralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEuropa = ( centralInstantAccelerationEuropa - centralInstantAccelerationJupiter ); //dependentVariables.begin( )->second.segment( 3, 3 ); //.transpose( );

//    std::cout << "initial acceleration Europa: " << centralInstantAccelerationEuropa.transpose( ) << "\n\n";
//    std::cout << "validation: " << centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 3, 3 ).transpose( ) << "\n\n";
//    std::cout << "test initial acceleration Europa: " << ( centralInstantAccelerationEuropa - centralInstantAccelerationJupiter ).transpose( ) << "\n\n";

//    // Compute acceleration partials.
//    Eigen::Matrix3d partialsAccelerationIoWrtPositionIo = Eigen::Matrix3d::Zero( );
//    partialsAccelerationIoWrtPositionIo += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialsAccelerationIoWrtPositionIo += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialsAccelerationIoWrtPositionIo += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialsAccelerationIoWrtPositionIo += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    Eigen::Matrix3d partialscentralInstantAccelerationEuropa = Eigen::Matrix3d::Zero( );
//    partialscentralInstantAccelerationEuropa += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialscentralInstantAccelerationEuropa += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialscentralInstantAccelerationEuropa += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    partialscentralInstantAccelerationEuropa += calculatePartialOfPointMassGravityFromRelativePositionWrtPositionOfAcceleratedBody(
//                relativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    // Compute acceleration partials.
    Eigen::Matrix3d partialsAccelerationIoWrtPositionIo;
    partialsAccelerationIoWrtPositionIo.block( 0, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 21, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 147, 3 ).transpose( );
    partialsAccelerationIoWrtPositionIo.block( 1, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 27, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 153, 3 ).transpose( );
    partialsAccelerationIoWrtPositionIo.block( 2, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 33, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 159, 3 ).transpose( );
//    partialsAccelerationIoWrtPositionIo = partialsAccelerationIoWrtPositionIo;
//    std::cout << "initial acceleration partials Io w.r.t. Io: " << partialsAccelerationIoWrtPositionIo << "\n\n";
//    std::cout << "same from dependent variables: " << partialsAccelerationIoWrtPositionIo << "\n\n";

    Eigen::Matrix3d partialsAccelerationIoWrtPositionEarth;
    partialsAccelerationIoWrtPositionEarth.block( 0, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 39, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 129, 3 ).transpose( );
    partialsAccelerationIoWrtPositionEarth.block( 1, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 45, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 135, 3 ).transpose( );
    partialsAccelerationIoWrtPositionEarth.block( 2, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 51, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 141, 3 ).transpose( );


    Eigen::Matrix3d partialsAccelerationEuropaWrtPositionEuropa;
    partialsAccelerationEuropaWrtPositionEuropa.block( 0, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 75, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 165, 3 ).transpose( );
    partialsAccelerationEuropaWrtPositionEuropa.block( 1, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 81, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 171, 3 ).transpose( );
    partialsAccelerationEuropaWrtPositionEuropa.block( 2, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 87, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 177, 3 ).transpose( );
//    partialscentralInstantAccelerationEuropa = partialsAccelerationEuropaWrtPositionEuropa;
//    std::cout << "initial acceleration partials Europa w.r.t. Europa: " << partialscentralInstantAccelerationEuropa << "\n\n";
//    std::cout << "same from dependent variables: " << partialsAccelerationEuropaWrtPositionEuropa << "\n\n";


    Eigen::Matrix3d partialsAccelerationEuropaWrtPositionEarth;
    partialsAccelerationEuropaWrtPositionEarth.block( 0, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 93, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 129, 3 ).transpose( );
    partialsAccelerationEuropaWrtPositionEarth.block( 1, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 99, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 135, 3 ).transpose( );
    partialsAccelerationEuropaWrtPositionEarth.block( 2, 0, 1, 3 ) = centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 105, 3 ).transpose( )
            - centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 141, 3 ).transpose( );

//    std::cout << "partialsAccelerationIoWrtPositionEarth: " << partialsAccelerationIoWrtPositionEarth << "\n\n";
//    std::cout << "partialsAccelerationEuropaWrtPositionEarth: " << partialsAccelerationEuropaWrtPositionEarth << "\n\n";


    // Compute initial acceleration Earth.
    Eigen::Vector3d relativePositionEarthSun = - sunCenteredState.segment( 0, 3 );
    Eigen::Vector3d relativePositionEarthJupiter = sunCenteredState.segment( 6, 3 ) - sunCenteredState.segment( 0, 3 );
    Eigen::Vector3d relativePositionEarthIo = sunCenteredState.segment( 12, 3 ) - sunCenteredState.segment( 0, 3 );
    Eigen::Vector3d relativePositionEarthEuropa = sunCenteredState.segment( 18, 3 ) - sunCenteredState.segment( 0, 3 );

    Eigen::Vector3d centralInstantAccelerationEarth = centralInstantDependentVariables.segment( 6, 3 ); //Eigen::Vector3d::Zero( );
//    centralInstantAccelerationEarth += calculatePointMassGravityAcceleration( relativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEarth += calculatePointMassGravityAcceleration( relativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEarth += calculatePointMassGravityAcceleration( relativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEarth += calculatePointMassGravityAcceleration( relativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    centralInstantAccelerationEarth = dependentVariables.begin( )->second.segment( 6, 3 ); //.transpose( );

//    std::cout << "initial acceleration Earth: " << centralInstantAccelerationEarth.transpose( ) << "\n\n";
//    std::cout << "validation: " << centralInstantDependentVariables/*dependentVariables.begin( )->second*/.segment( 6, 3 ).transpose( ) << "\n\n";

    Eigen::Vector3d partialOfSecondTimeDerivativeRightAscensionWrtIoPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth,
                                                                                  partialsAccelerationIoWrtPositionIo, true );

    Eigen::Vector3d partialOfSecondTimeDerivativeDeclinationWrtIoPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
                                                                               ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                                                                               centralInstantAccelerationIo - centralInstantAccelerationEarth,
                                                                               partialsAccelerationIoWrtPositionIo, true );

    double secondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                centralInstantAccelerationIo - centralInstantAccelerationEarth );

    double secondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                centralInstantAccelerationIo - centralInstantAccelerationEarth );

    Eigen::Vector2d instrumentalFrameRelativeAcceleration = mutualApproximationScaling->getInstrumentalFrameRelativeAcceleration( );

    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtPositionIo =
            mutualApproximationScaling->getPartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition( );
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEuropa =
            mutualApproximationScaling->getPartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition( );
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEarth =
            mutualApproximationScaling->getPartialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition( );

    Eigen::Vector3d partialOfSecondTimeDerivativeXwrtPositionIo = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionIo.block( 0, 0, 1, 3 ).transpose( );
//            computepartialOfSecondTimeDerivativeXwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                partialsAccelerationIoWrtPositionIo, partialsAccelerationEuropaWrtPositionEuropa, true );

    Eigen::Vector3d partialOfSecondTimeDerivativeXwrtPositionEarth = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEarth.block( 0, 0, 1, 3 ).transpose( ); //- 1.0 * (
//            computepartialOfSecondTimeDerivativeXwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true ) );

    Eigen::Vector3d partialOfSecondTimeDerivativeXwrtPositionEuropa = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEuropa.block( 0, 0, 1, 3 ).transpose( );
//            computepartialOfSecondTimeDerivativeXwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           partialsAccelerationIoWrtPositionIo, partialsAccelerationEuropaWrtPositionEuropa, false );

//    partialOfSecondTimeDerivativeXwrtPositionEarth -=
//            computepartialOfSecondTimeDerivativeXwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );

    double initialSecondPartialXwrtTime = instrumentalFrameRelativeAcceleration[ 0 ]; //computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                        ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                        centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                        rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                                                        rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

    Eigen::Vector3d partialOfSecondTimeDerivativeYwrtPositionIo = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionIo.block( 1, 0, 1, 3 ).transpose( );
//            computepartialOfSecondTimeDerivativeYwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

    Eigen::Vector3d partialOfSecondTimeDerivativeYwrtPositionEarth = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEarth.block( 1, 0, 1, 3 ).transpose( ); //- 1.0 * (
//            computepartialOfSecondTimeDerivativeYwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true ) );

    Eigen::Vector3d partialOfSecondTimeDerivativeYwrtPositionIoForEuropa = partialsOfInstrumentalFrameRelativeAccelerationWrtPositionEuropa.block( 1, 0, 1, 3 ).transpose( );
//            computepartialOfSecondTimeDerivativeYwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false );

//    partialOfSecondTimeDerivativeYwrtPositionEarth -=
//            computepartialOfSecondTimeDerivativeYwrtPositionIo( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                           - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );


    double initialSecondPartialYwrtTime = instrumentalFrameRelativeAcceleration[ 1 ]; // computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                        ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                        centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );



    /// x - axis
    Eigen::Vector3d variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << variationCentralInstantPositionIoEarth[ 0 ], 0.0, 0.0 ).finished( );
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +  variationCentralInstantPositionIo;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    Eigen::Vector3d modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
    Eigen::Vector3d modifiedRelativePositionIoEarth = relativePositionIoEarth - variationCentralInstantPositionIo;
    Eigen::Vector3d modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
    Eigen::Vector3d modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;


//    // Define propagator settings.
//    Eigen::VectorXd modifiedInitialState = systemInitialState;
//    Eigen::VectorXd variationStateCentralInstant = Eigen::VectorXd::Zero( 24 );
//    variationStateCentralInstant[ 12 ] = variationCentralInstantPositionIoEarth[ 0 ];

//    Eigen::VectorXd variationInitialState = stateTransitionInterpolator->interpolate( estimatedCentralInstant ).inverse( ) * variationStateCentralInstant;
//    modifiedInitialState += variationInitialState;

//    propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
//                centralBodies, accelerationModelMap, bodiesToPropagate, modifiedInitialState, std::make_shared< propagators::PropagationTimeTerminationSettings > ( simulationEndEpoch, true ),
//                propagators::cowell, dependentVariablesToSave );

//    // Create simulation object and propagate dynamics.
//    propagators::SingleArcVariationalEquationsSolver< >
//            modifiedVariationalEquationsSimulator( bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
//                                                   std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
//                                                   dependentVariablesToSave );


    // Define propagator settings.
    Eigen::VectorXd modifiedInitialState = systemInitialState;
    Eigen::VectorXd variationStateCentralInstant = Eigen::VectorXd::Zero( 24 );
    variationStateCentralInstant[ 12 ] = variationCentralInstantPositionIo[ 0 ];

//    std::cout << "variation state central state: " << variationStateCentralInstant.transpose( ) << "\n\n";

    Eigen::VectorXd variationInitialState = stateTransitionInterpolator->interpolate( estimatedCentralInstant ).inverse( ) * variationStateCentralInstant;

    std::cout << "test inverted state transition matrix: " <<
                 ( variationStateCentralInstant - stateTransitionInterpolator->interpolate( estimatedCentralInstant ) * variationInitialState ).transpose( ) << "\n\n";

    modifiedInitialState += variationInitialState;
//    std::cout << "variation initial state: " << variationInitialState.transpose( ) << "\n\n";

//    modifiedInitialState.segment( 12, 3 ) += stateTransitionInterpolator->interpolate( estimatedCentralInstant ).inverse( ) * variationCentralInstantPositionIo;
    propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, modifiedInitialState, std::make_shared< propagators::PropagationTimeTerminationSettings > ( simulationEndEpoch, true ),
                propagators::cowell, dependentVariablesToSave );

    // Create simulation object and propagate dynamics.
    propagators::SingleArcVariationalEquationsSolver< > modifiedVariationalEquationsSimulator( bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true, dependentVariablesToSave );

    std::map< double, Eigen::VectorXd > modifiedStateHistory = modifiedVariationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > modifiedDependentVariablesHistory = modifiedVariationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > modifiedStateHistoryInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( modifiedStateHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( modifiedStateHistory ), 4 );


    // Update apparent distance observation model.
    apparentDistanceModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, observableSettings, bodyMap );

    apparentDistanceHistory.clear( );
    for ( std::map< double, Eigen::VectorXd >::iterator itr = modifiedStateHistory.begin( ) ; itr!= modifiedStateHistory.end( ) ; itr++ )
    {
        Eigen::Vector3d positionIoEarth = modifiedDependentVariablesHistory[ itr->first ].segment( 9, 3 );
        Eigen::Vector3d positionEuropaEarth = modifiedDependentVariablesHistory[ itr->first ].segment( 15, 3 );
        std::pair< double, double > angularPositionIoWrtEarth = computeRightAscensionAndDeclination( positionIoEarth );
        std::pair< double, double > angularPositionEuropaWrtEarth = computeRightAscensionAndDeclination( positionEuropaEarth );
        apparentDistanceHistory[ itr->first ] = (
                    Eigen::Vector1d( ) << 3600.0 * 180.0 / mathematical_constants::PI * /*apparentDistanceModel->computeObservations( itr->first, observation_models::receiver )*/
                    std::sqrt( ( angularPositionEuropaWrtEarth.first - angularPositionIoWrtEarth.first )
                               * std::cos( ( angularPositionEuropaWrtEarth.second + angularPositionIoWrtEarth.second ) / 2.0 )
                               * ( angularPositionEuropaWrtEarth.first - angularPositionIoWrtEarth.first )
                               * std::cos( ( angularPositionEuropaWrtEarth.second + angularPositionIoWrtEarth.second ) / 2.0 )
                               + ( angularPositionEuropaWrtEarth.second - angularPositionIoWrtEarth.second )
                               * ( angularPositionEuropaWrtEarth.second - angularPositionIoWrtEarth.second ) ) ).finished( );
    }


    input_output::writeDataMapToTextFile( apparentDistanceHistory,
                                          "apparentDistanceHistory_variationXio.dat",
                                          "C:/Users/chamb/Documents/PhD/testMutualApproximationPartials/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );



//    double newEstimatedCentralInstant = computeCentralInstantFromDependentVariables( simulationStartEpoch, 60.0,
//                                                                                     modifiedDependentVariablesHistory );
//    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
//    std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";
//    std::cout << "numerical partial central instant w.r.t. first transmitter position: " << ( newEstimatedCentralInstant - estimatedCentralInstant ) / ( variationCentralInstantPositionIo[ 0 ] ) << "\n\n";
//    std::cout << "variation along x-axis for first transmitter position: " <<  variationCentralInstantPositionIo[ 0 ] << "\n\n";

//    std::cout << "modified state at central instant: " << ( modifiedStateHistoryInterpolator->interpolate( estimatedCentralInstant ) - variationStateCentralInstant ).transpose( ) << "\n\n";
//    std::cout << "initial state at central instant: " << stateHistoryInterpolator->interpolate( estimatedCentralInstant ).transpose( ) << "\n\n";
//    std::cout << "difference: " << ( modifiedStateHistoryInterpolator->interpolate( estimatedCentralInstant ) - variationStateCentralInstant
//                                     - stateHistoryInterpolator->interpolate( estimatedCentralInstant ) ).transpose( ) << "\n\n";


//    // Retrieve dependent variables history.
//    std::map< double, Eigen::VectorXd > modifiedDependentVariablesHistory = modifiedVariationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );
//    // Create dependent variables interpolator.
//    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > modifiedDependentVariablesInterpolator
//            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
//                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( modifiedDependentVariablesHistory ),
//                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( modifiedDependentVariablesHistory ), 4 );

//    std::shared_ptr< propagators::SingleArcDependentVariablesInterface > modifiedDependentVariablesInterface = std::make_shared< propagators::SingleArcDependentVariablesInterface >(
//                modifiedDependentVariablesInterpolator, dependentVariablesToSave );

//    modifiedMutualApproximationScaling =
//            std::make_shared< observation_partials::MutualApproximationScaling >( modifiedDependentVariablesInterface );


//    modifiedLinkEndStates.clear( );
//    modifiedCentralInstantIo = centralInstantIoState;
//    modifiedCentralInstantIo[ 0 ] += variationCentralInstantPositionIoEarth[ 0 ];
//    modifiedLinkEndStates.push_back( modifiedCentralInstantIo );
//    modifiedLinkEndStates.push_back( centralInstantEuropaState );
//    modifiedLinkEndStates.push_back( centralInstantEarthState );
//    modifiedMutualApproximationScaling->update( modifiedLinkEndStates, linkEndsTimes, observation_models::receiver, linkEnds, modifiedObservation );


    // TO BE MODIFIED AND RETRIEVED FROM DEPENDENT VARIABLES
    Eigen::Vector3d modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun,
                                                                                   bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth,
                                                                                   bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter,
                                                                                   bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa,
                                                                                   bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    Eigen::Vector3d modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
    Eigen::Vector3d modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun,
                                                                                        bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo,
                                                                                        bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa,
                                                                                        bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

    double modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    double variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of second partial right ascension wrt time - x: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeRightAscensionWrtIoPosition[ 0 ] << "\n\n";


    double modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    double variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of second partial declination wrt time - x: " << variationSecondPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeDeclinationWrtIoPosition[ 0 ] << "\n\n";

//    std::cout << "partial acceleration wrt variation in x: " << ( ( modifiedCentralInstantAccelerationIo - centralInstantAccelerationIo )
//                                                                  / variationCentralInstantPositionIoEarth[ 0 ] ).transpose( ) << "\n\n";
//    std::cout << "validation: " << partialsAccelerationIoWrtPositionIo.block( 0, 0, 3, 1 ).transpose( ) << "\n\n";


    // CANNOT USE INSTRUMENTAL FRAM ACCELERATION FROM MODIFIED MUTUAL APPROXIMATION SCALING
    // BECAUSE THE DEPENDENT VARIABLES ARE NOT UPDATED -> ISSUE WITH THE ACCELERATION VALUES
    Eigen::Vector2d modifiedInstrumentalFrameRelativeAcceleration = modifiedMutualApproximationScaling->getInstrumentalFrameRelativeAcceleration( );
    double modifiedSecondTimeDerivativeOfX = /*modifiedInstrumentalFrameRelativeAcceleration[ 0 ];*/ computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                         ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                         modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                         modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                         rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    double variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionIoEarth[ 0 ];

    double modifiedSecondTimeDerivativeOfY = /*modifiedInstrumentalFrameRelativeAcceleration[ 1 ];*/ computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                         ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                         modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    double variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionIoEarth[ 0 ];

    std::cout << "numerical position partial of second partial X wrt time - x: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeXwrtPositionIo[ 0 ] << "\n\n";
    std::cout << "numerical position partial of second partial Y wrt time - x: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeYwrtPositionIo[ 0 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - x (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionIo[ 0 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - x (model): " << ( -2.137060741744073e-20 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - x (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIo[ 0 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - x (model): " << ( -2.595182334963883e-21 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";

    double numericalPartialCentralInstantWrtIoPosition = - 0.349823200512844 / variationCentralInstantPositionIoEarth[ 0 ];
    std::cout << "numerical partial of central instant w.r.t. Io position - x: " << numericalPartialCentralInstantWrtIoPosition << "\n\n";


    ///// W.R.T. SECOND TRANSMITTER POSITION, I.E. W.R.T. Europa position
    Eigen::Vector3d variationInitialVectorEuropa = ( Eigen::Vector3d( ) << variationCentralInstantPositionEuropaEarth[ 0 ], 0.0, 0.0 ).finished( );
    Eigen::Vector3d modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
    std::pair< double, double > modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
    Eigen::Vector3d modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
    Eigen::Vector3d modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
    Eigen::Vector3d modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
    Eigen::Vector3d modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

    Eigen::Vector3d modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    Eigen::Vector3d modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 0 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 0 ];

    std::cout << "Europa: numerical position partial of second partial X wrt time - x: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeXwrtPositionEuropa[ 0 ] << "\n\n";
    std::cout << "Europa: numerical position partial of second partial Y wrt time - x: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 0 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - x (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEuropa[ 0 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - x (model): " << ( 3.978408922438975e-21 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - x (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 0 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - x (model): " << ( 5.212567519655293e-22 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";



    ///// W.R.T. RECEIVER POSITION, I.E. W.R.T. Earth position
    Eigen::Vector3d variationInitialVectorEarth = ( Eigen::Vector3d( ) << variationCentralInstantPositionEuropaEarth[ 0 ], 0.0, 0.0 ).finished( );

    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo; // - variationInitialVectorEuropa;

    Eigen::Vector3d modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
    Eigen::Vector3d modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
    Eigen::Vector3d modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
    Eigen::Vector3d modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

    Eigen::Vector3d modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );


//    //-----------------
//    // TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE

//    Eigen::Vector3d positionPartialOfSecondTimeDerivativeOfRightAscensionWrtReceiverState = -
//            computepartialOfSecondTimeDerivativeRightAscensionWrtIoPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                       ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                        centralInstantAccelerationIo - centralInstantAccelerationEarth, - partialsAccelerationIoWrtPositionEarth );

//    modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 0 ];

//    std::cout << "--- Test position partial of second time derivative of right ascension w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of right ascension wrt RECEIVER STATE - x: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfRightAscensionWrtReceiverState[ 0 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    Eigen::Vector3d positionPartialOfSecondTimeDerivativeOfDeclinationWrtReceiverState = -
//            computepartialOfSecondTimeDerivativeDeclinationWrtIoPosition( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                     ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                     centralInstantAccelerationIo - centralInstantAccelerationEarth, - partialsAccelerationIoWrtPositionEarth );

//    modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 0 ];

//    std::cout << "--- Test position partial of second time derivative of declination w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of declination wrt RECEIVER STATE - x: " << variationSecondPartialDeclinationWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfDeclinationWrtReceiverState[ 0 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    // END TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE
//    //-----------------

    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 0 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 0 ];

    std::cout << "Earth: numerical position partial of second partial X wrt time - x: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeXwrtPositionEarth[ 0 ] << "\n\n";
    std::cout << "Earth: numerical position partial of second partial Y wrt time - x: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeYwrtPositionEarth[ 0 ] << "\n\n";


//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - x (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEarth[ 0 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - x (model): " << ( 1.902079477666416e-024 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - x (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionEarth[ 0 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - x (model): " << ( -1.799037379919735e-023 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";




    /// y - axis
    variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << 0.0,  variationCentralInstantPositionIoEarth[ 1 ], 0.0 ).finished( );
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +  variationCentralInstantPositionIo;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationCentralInstantPositionIo;
    modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
    modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;

    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

    modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of second partial right ascension wrt time - y: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeRightAscensionWrtIoPosition[ 1 ] << "\n\n";


    modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of second partial declination wrt time - y: " << variationSecondPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeDeclinationWrtIoPosition[ 1 ] << "\n\n";


    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                         ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                         modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                         modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                         rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionIoEarth[ 1 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionIoEarth[ 1 ];

    std::cout << "numerical position partial of second partial X wrt time - y: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeXwrtPositionIo[ 1 ] << "\n\n";
    std::cout << "numerical position partial of second partial Y wrt time - y: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeYwrtPositionIo[ 1 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - y (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionIo[ 1 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - y (model): " << ( -6.937614422123296e-021 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - y (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIo[ 1 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - y (model): " << ( -5.453963675201478e-021 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";


    ///// W.R.T. SECOND TRANSMITTER POSITION, I.E. W.R.T. Europa position
    variationInitialVectorEuropa = ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionEuropaEarth[ 1 ], 0.0 ).finished( );
    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
    modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 1 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 1 ];

    std::cout << "Europa: numerical position partial of second partial X wrt time - y: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeXwrtPositionEuropa[ 1 ] << "\n\n";
    std::cout << "Europa: numerical position partial of second partial Y wrt time - y: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 1 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - y (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEuropa[ 1 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - y (model): " << ( 3.321504421999323e-021 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - y (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 1 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - y (model): " << ( 1.548071925631198e-021 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";


    ///// W.R.T. RECEIVER POSITION, I.E. W.R.T. Earth position
    variationInitialVectorEarth = ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionEuropaEarth[ 1 ], 0.0 ).finished( );

    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEuropa;
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEuropa;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo; // - variationInitialVectorEuropa;

    modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
    modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
    modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
    modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

    modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );


//    //-----------------
//    // TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE

//    modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 1 ];

//    std::cout << "--- Test position partial of second time derivative of right ascension w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of right ascension wrt RECEIVER STATE - y: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfRightAscensionWrtReceiverState[ 1 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 1 ];

//    std::cout << "--- Test position partial of second time derivative of declination w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of declination wrt RECEIVER STATE - y: " << variationSecondPartialDeclinationWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfDeclinationWrtReceiverState[ 1 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    // END TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE
//    //-----------------


    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 1 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 1 ];

    std::cout << "Earth: numerical position partial of second partial X wrt time - y: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeXwrtPositionEarth[ 1 ] << "\n\n";
    std::cout << "Earth: numerical position partial of second partial Y wrt time - y: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeYwrtPositionEarth[ 1 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - y (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEarth[ 1 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - y (model): " << ( -4.075550828425433e-023 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - y (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionEarth[ 1 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - y (model): " << ( -4.911865021105991e-024 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";


    /// z- axis
    variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionIoEarth[ 2 ] ).finished( );
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) + variationCentralInstantPositionIo;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
    modifiedRelativePositionIoEarth = relativePositionIoEarth - variationCentralInstantPositionIo;
    modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
    modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;

    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

    modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
            / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of second partial right ascension wrt time - z: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeRightAscensionWrtIoPosition[ 2 ] << "\n\n";


    modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth );

    variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
            / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of second partial declination wrt time - z: " << variationSecondPartialDeclinationWrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeDeclinationWrtIoPosition[ 2 ] << "\n\n";


    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                  rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionIoEarth[ 2 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionIoEarth[ 2 ];

    std::cout << "numerical position partial of second partial X wrt time - z: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeXwrtPositionIo[ 2 ] << "\n\n";
    std::cout << "numerical position partial of second partial Y wrt time - z: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "validation: " << partialOfSecondTimeDerivativeYwrtPositionIo[ 2 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - z (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionIo[ 2 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of X - z (model): " << ( -4.988650059491025e-021 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - z (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIo[ 2 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Io: relative diff for position partials of second time derivative of Y - z (model): " << ( 1.016433129881544e-020 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";


    ///// W.R.T. SECOND TRANSMITTER POSITION, I.E. W.R.T. Europa position
    variationInitialVectorEuropa = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionEuropaEarth[ 2 ] ).finished( );
    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
    modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
    modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

    modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 2 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 2 ];

    std::cout << "Europa: numerical position partial of second partial X wrt time - z: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeXwrtPositionEuropa[ 2 ] << "\n\n";
    std::cout << "Europa: numerical position partial of second partial Y wrt time - z: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Europa: validation: " << partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 2 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - z (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEuropa[ 2 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of X - z (model): " << ( 1.941445730071036e-021 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - z (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionIoForEuropa[ 2 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Europa: relative diff for position partials of second time derivative of Y - z (model): " << ( -2.48101951029586e-021 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";

    ///// W.R.T. RECEIVER POSITION, I.E. W.R.T. Earth position
    variationInitialVectorEarth = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionEuropaEarth[ 2 ] ).finished( );

    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEuropa;
    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEuropa;
    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter; // - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo; // - variationInitialVectorEuropa;

    modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
    modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
    modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
    modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

    modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    //-----------------
//    // TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE

//    modifiedSecondPartialRightAscensionWrtTime = computeSecondPartialRightAscensionWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialRightAscensionWrtTime = ( modifiedSecondPartialRightAscensionWrtTime - secondPartialRightAscensionWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 2 ];

//    std::cout << "--- Test position partial of second time derivative of right ascension w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of right ascension wrt RECEIVER STATE - z: " << variationSecondPartialRightAscensionWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfRightAscensionWrtReceiverState[ 2 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    modifiedSecondPartialDeclinationWrtTime = computeSecondPartialDeclinationWrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth );

//    variationSecondPartialDeclinationWrtTime = ( modifiedSecondPartialDeclinationWrtTime - secondPartialDeclinationWrtTime )
//            / variationCentralInstantPositionEuropaEarth[ 2 ];

//    std::cout << "--- Test position partial of second time derivative of declination w.r.t. RECEIVER state ---" << "\n\n";
//    std::cout << "numerical position partial of second time derivative of declination wrt RECEIVER STATE - z: " << variationSecondPartialDeclinationWrtTime << "\n\n";
//    std::cout << "validation: " << positionPartialOfSecondTimeDerivativeOfDeclinationWrtReceiverState[ 2 ] << "\n\n";
//    std::cout << "---" << "\n\n";

//    // END TEST POSITION PARTIAL OF SECOND TIME DERIVATIVE RIGHT ASCENSION/DECLINATION W.R.T. RECEIVER STATE
//    //-----------------

    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );
    variationSecondPartialXwrtTime = ( modifiedSecondTimeDerivativeOfX - initialSecondPartialXwrtTime ) / variationCentralInstantPositionEuropaEarth[ 2 ];

    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
                                                                  centralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, centralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );
    variationSecondPartialYwrtTime = ( modifiedSecondTimeDerivativeOfY - initialSecondPartialYwrtTime ) / variationCentralInstantPositionEuropaEarth[ 2 ];

    std::cout << "Earth: numerical position partial of second partial X wrt time - z: " << variationSecondPartialXwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeXwrtPositionEarth[ 2 ] << "\n\n";
    std::cout << "Earth: numerical position partial of second partial Y wrt time - z: " << variationSecondPartialYwrtTime << "\n\n";
    std::cout << "Earth: validation: " << partialOfSecondTimeDerivativeYwrtPositionEarth[ 2 ] << "\n\n";

//    std::cout << "------------------------------" << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - z (functions): " << ( partialOfSecondTimeDerivativeXwrtPositionEarth[ 2 ] -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of X - z (model): " << ( -1.735991754417374e-023 -  variationSecondPartialXwrtTime ) / variationSecondPartialXwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - z (functions): " << ( partialOfSecondTimeDerivativeYwrtPositionEarth[ 2 ] -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "Earth: relative diff for position partials of second time derivative of Y - z (model): " << ( -9.915275817070896e-024 -  variationSecondPartialYwrtTime ) / variationSecondPartialYwrtTime << "\n\n";
//    std::cout << "------------------------------" << "\n\n";


//    std::cout << "total acceleration partial of Earth w.r.t. state of Earth: " << dependentVariablesHistory.rbegin( )->second.segment( 93, 18 ).transpose( ) << "\n\n";
//    std::cout << "total acceleration partial of Earth w.r.t. state of Io: " << dependentVariablesHistory.rbegin( )->second.segment( 111, 18 ).transpose( ) << "\n\n";
//    std::cout << "total acceleration partial of Earth w.r.t. state of Europa: " << dependentVariablesHistory.rbegin( )->second.segment( 129, 18 ).transpose( ) << "\n\n";



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////        TEST POSITION PARTIALS OF INTERMEDIATE VARIABLES USED TO DERIVE THE CENTRAL INSTANT         ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    std::cout << "----------------------------------------------------------------------------------------------------------------------" << "\n\n";
//    std::cout << "--------------- TEST POSITION PARTIALS OF INTERMEDIATE VARIABLES USED TO DERIVE THE CENTRAL INSTANT ------------------" << "\n\n";
//    std::cout << "----------------------------------------------------------------------------------------------------------------------" << "\n\n";

//    // Compute coefficients of the cubic polynomial for the central instant t0.
//    Eigen::Vector4d cubicPolynomialCoefficients =
//            ( Eigen::Vector4d( ) <<

//              initialSecondPartialXwrtTime * initialSecondPartialXwrtTime
//            + initialSecondPartialYwrtTime * initialSecondPartialYwrtTime,

//            3.0 * ( initialSecondPartialXwrtTime * firstTimePartialOfX
//            + initialSecondPartialYwrtTime * firstTimePartialOfY ),

//            2.0 * ( centralInstantValueX * initialSecondPartialXwrtTime
//            + centralInstantValueY * initialSecondPartialYwrtTime
//            + firstTimePartialOfX * firstTimePartialOfX
//            + firstTimePartialOfY * firstTimePartialOfY ),

//            2.0 * ( centralInstantValueX * firstTimePartialOfX
//            + centralInstantValueY * firstTimePartialOfY ) ).finished( );

//    std::cout << "cubicPolynomialCoefficients (function): " << cubicPolynomialCoefficients.transpose( ) << "\n\n";

//    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
//    Eigen::Vector3d depressedCubicPolynomialCoefficients =
//            ( Eigen::Vector3d( ) <<
//              cubicPolynomialCoefficients[ 1 ] / cubicPolynomialCoefficients[ 0 ],
//            cubicPolynomialCoefficients[ 2 ] / cubicPolynomialCoefficients[ 0 ],
//            cubicPolynomialCoefficients[ 3 ] / cubicPolynomialCoefficients[ 0 ] ).finished( );

//    double Q = computeConstantQ( depressedCubicPolynomialCoefficients[ 1 ], depressedCubicPolynomialCoefficients[ 0 ] );
//    double R = computeConstantR( depressedCubicPolynomialCoefficients[ 2 ], depressedCubicPolynomialCoefficients[ 1 ], depressedCubicPolynomialCoefficients[ 0 ] );

//    std::cout << "intermediate variable Q (function): " << Q << "\n\n";
//    std::cout << "intermediate variable R (function): " << R << "\n\n";


//    double angleTheta = computeAngleThetaSolutionCubicEquation(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );



//    double beta = Q * Q * Q + R * R;

//    double S = 0.0;
//    if ( ( R + sqrt( beta ) ) >= 0 )
//    {
//        S = std::pow( R + sqrt( beta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        S = - std::pow( std::fabs( R + std::sqrt( beta ) ), 1.0 / 3.0 );
//    }

//    double T = 0.0;
//    if ( ( R - sqrt( beta ) ) >= 0 )
//    {
//        T = std::pow( R - sqrt( beta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        T = - std::pow( std::fabs( R - std::sqrt( beta ) ), 1.0 / 3.0 );
//    }



//    Eigen::Vector3d numericalPartialsQwrtIoPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsQwrtEuropaPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsQwrtEarthPosition = Eigen::Vector3d::Zero( );

//    Eigen::Vector3d numericalPartialsRwrtIoPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsRwrtEuropaPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsRwrtEarthPosition = Eigen::Vector3d::Zero( );

//    Eigen::Vector3d numericalPartialsThetawrtIoPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsThetawrtEuropaPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsThetawrtEarthPosition = Eigen::Vector3d::Zero( );

//    Eigen::Vector3d numericalPartialsSwrtIoPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsSwrtEuropaPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsSwrtEarthPosition = Eigen::Vector3d::Zero( );

//    Eigen::Vector3d numericalPartialsTwrtIoPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsTwrtEuropaPosition = Eigen::Vector3d::Zero( );
//    Eigen::Vector3d numericalPartialsTwrtEarthPosition = Eigen::Vector3d::Zero( );

//    /// x - axis
//    /// w.r.t. position of first transmitter
//    variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << variationCentralInstantPositionIoEarth[ 0 ], 0.0, 0.0 ).finished( );
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +  variationCentralInstantPositionIo;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEarth = relativePositionIoEarth - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );


//    Eigen::Vector4d modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    Eigen::Vector3d modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    double modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    double modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtIoPosition[ 0 ] = ( modifiedQ - Q ) / variationCentralInstantPositionIoEarth[ 0 ];
//    numericalPartialsRwrtIoPosition[ 0 ] = ( modifiedR - R ) / variationCentralInstantPositionIoEarth[ 0 ];


//    double modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtIoPosition[ 0 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionIoEarth[ 0 ];


//    double modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    double modifiedS = 0.0;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    double modifiedT = 0.0;
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtIoPosition[ 0 ] = ( modifiedS - S ) / variationCentralInstantPositionIoEarth[ 0 ];
//    numericalPartialsTwrtIoPosition[ 0 ] = ( modifiedT - T ) / variationCentralInstantPositionIoEarth[ 0 ];


//    // Define propagator settings.
//    Eigen::VectorXd modifiedInitialState = systemInitialState;
//    Eigen::VectorXd variationStateCentralInstant = Eigen::VectorXd::Zero( 24 );
//    variationStateCentralInstant[ 12 ] = 0.0001 * currentPropagatedState[ 12 ]; //variationCentralInstantPositionIo;

////    std::cout << "variation state central state: " << variationStateCentralInstant.transpose( ) << "\n\n";

//    Eigen::VectorXd variationInitialState = stateTransitionInterpolator->interpolate( estimatedCentralInstant ).inverse( ) * variationStateCentralInstant;

//    std::cout << "test inverted state transition matrix: " <<
//                 ( variationStateCentralInstant - stateTransitionInterpolator->interpolate( estimatedCentralInstant ) * variationInitialState ).transpose( ) << "\n\n";

//    modifiedInitialState += variationInitialState;
////    std::cout << "variation initial state: " << variationInitialState.transpose( ) << "\n\n";

////    modifiedInitialState.segment( 12, 3 ) += stateTransitionInterpolator->interpolate( estimatedCentralInstant ).inverse( ) * variationCentralInstantPositionIo;
//    propagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > > (
//                centralBodies, accelerationModelMap, bodiesToPropagate, modifiedInitialState, std::make_shared< propagators::PropagationTimeTerminationSettings > ( simulationEndEpoch, true ),
//                propagators::cowell, dependentVariablesToSave );

//    // Create simulation object and propagate dynamics.
//    SingleArcVariationalEquationsSolver< > modifiedVariationalEquationsSimulator( bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true, dependentVariablesToSave );

//    std::map< double, Eigen::VectorXd > modifiedStateHistory = modifiedVariationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > modifiedDependentVariablesHistory = modifiedVariationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );
//    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > modifiedStateHistoryInterpolator
//            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
//                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( modifiedStateHistory ),
//                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( modifiedStateHistory ), 4 );

//    double newEstimatedCentralInstant = computeCentralInstantFromDependentVariables( simulationStartEpoch, 60.0,
//                                                                                     modifiedDependentVariablesHistory );
//    std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";
//    std::cout << "numerical partial central instant w.r.t. first transmitter position: " << ( newEstimatedCentralInstant - estimatedCentralInstant ) / ( 0.0001 * currentPropagatedState[ 12 ] ) << "\n\n";
//    std::cout << "variation along x-axis for first transmitter position: " <<  0.0001 * currentPropagatedState[ 12 ] << "\n\n";

//    std::cout << "modified state at central instant: " << ( modifiedStateHistoryInterpolator->interpolate( estimatedCentralInstant ) - variationStateCentralInstant ).transpose( ) << "\n\n";
//    std::cout << "initial state at central instant: " << stateHistoryInterpolator->interpolate( estimatedCentralInstant ).transpose( ) << "\n\n";
//    std::cout << "difference: " << ( modifiedStateHistoryInterpolator->interpolate( estimatedCentralInstant ) - variationStateCentralInstant
//                                     - stateHistoryInterpolator->interpolate( estimatedCentralInstant ) ).transpose( ) << "\n\n";


////    // Create observation model.
////    observationModel = observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
////                linkEnds, observableSettings, bodyMap );

////    Eigen::Vector1d modifiedCentralInstant = observationModel->computeObservations( estimatedCentralInstant, observation_models::receiver );
////    std::cout << "modified central instant: " << modifiedCentralInstant[ 0 ] << "\n\n";
////    std::cout << "numerical partial central instant w.r.t. first transmitter position: " << ( modifiedCentralInstant[ 0 ] - estimatedCentralInstant ) / ( 0.0001 * variationCentralInstantPositionIo[ 0 ] ) << "\n\n";



//    /// w.r.t. position of second transmitter
//    variationInitialVectorEuropa = ( Eigen::Vector3d( ) << variationCentralInstantPositionEuropaEarth[ 0 ], 0.0, 0.0 ).finished( );
//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEuropaPosition[ 0 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 0 ];
//    numericalPartialsRwrtEuropaPosition[ 0 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 0 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEuropaPosition[ 0 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 0 ];

//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEuropaPosition[ 0 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 0 ];
//    numericalPartialsTwrtEuropaPosition[ 0 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 0 ];



//    /// w.r.t. position of receiver
//    variationInitialVectorEarth = ( Eigen::Vector3d( ) << variationCentralInstantPositionEuropaEarth[ 0 ], 0.0, 0.0 ).finished( );

//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
//    modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
//    modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
//    modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

//    modifiedRelativePositionIoEarth = relativePositionIoEarth + variationInitialVectorEarth;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth + variationInitialVectorEarth;

//    modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= centralInstantAccelerationJupiter;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= centralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEarthPosition[ 0 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 0 ];
//    numericalPartialsRwrtEarthPosition[ 0 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 0 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEarthPosition[ 0 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 0 ];

//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEarthPosition[ 0 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 0 ];
//    numericalPartialsTwrtEarthPosition[ 0 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 0 ];



//    /// y - axis
//    /// w.r.t. position of first transmitter
//    variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionIoEarth[ 1 ], 0.0 ).finished( );
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +  variationCentralInstantPositionIo;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEarth = relativePositionIoEarth - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );


//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtIoPosition[ 1 ] = ( modifiedQ - Q ) / variationCentralInstantPositionIoEarth[ 1 ];
//    numericalPartialsRwrtIoPosition[ 1 ] = ( modifiedR - R ) / variationCentralInstantPositionIoEarth[ 1 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtIoPosition[ 1 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionIoEarth[ 1 ];


//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtIoPosition[ 1 ] = ( modifiedS - S ) / variationCentralInstantPositionIoEarth[ 1 ];
//    numericalPartialsTwrtIoPosition[ 1 ] = ( modifiedT - T ) / variationCentralInstantPositionIoEarth[ 1 ];



//    /// w.r.t. position of second transmitter
//    variationInitialVectorEuropa = ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionEuropaEarth[ 1 ], 0.0 ).finished( );
//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEuropaPosition[ 1 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 1 ];
//    numericalPartialsRwrtEuropaPosition[ 1 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 1 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEuropaPosition[ 1 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 1 ];


//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEuropaPosition[ 1 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 1 ];
//    numericalPartialsTwrtEuropaPosition[ 1 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 1 ];



//    /// w.r.t. position of receiver
//    variationInitialVectorEarth = ( Eigen::Vector3d( ) << 0.0, variationCentralInstantPositionEuropaEarth[ 1 ], 0.0 ).finished( );

//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
//    modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
//    modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
//    modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

//    modifiedRelativePositionIoEarth = relativePositionIoEarth + variationInitialVectorEarth;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth + variationInitialVectorEarth;

//    modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= centralInstantAccelerationJupiter;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= centralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEarthPosition[ 1 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 1 ];
//    numericalPartialsRwrtEarthPosition[ 1 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 1 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEarthPosition[ 1 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 1 ];


//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEarthPosition[ 1 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 1 ];
//    numericalPartialsTwrtEarthPosition[ 1 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 1 ];




//    /// z - axis
//    /// w.r.t. position of first transmitter
//    variationCentralInstantPositionIo = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionIoEarth[ 2 ] ).finished( );
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) +  variationCentralInstantPositionIo;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedRelativePositionIoSun = relativePositionIoSun - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEarth = relativePositionIoEarth - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoJupiter = relativePositionIoJupiter - variationCentralInstantPositionIo;
//    modifiedRelativePositionIoEuropa = relativePositionIoEuropa - variationCentralInstantPositionIo;

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedRelativePositionJupiterIo = relativePositionJupiterIo + variationCentralInstantPositionIo;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedRelativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  rightAscensionAndDeclinationEuropa.first, rightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth );


//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtIoPosition[ 2 ] = ( modifiedQ - Q ) / variationCentralInstantPositionIoEarth[ 2 ];
//    numericalPartialsRwrtIoPosition[ 2 ] = ( modifiedR - R ) / variationCentralInstantPositionIoEarth[ 2 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtIoPosition[ 2 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionIoEarth[ 2 ];

//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtIoPosition[ 2 ] = ( modifiedS - S ) / variationCentralInstantPositionIoEarth[ 2 ];
//    numericalPartialsTwrtIoPosition[ 2 ] = ( modifiedT - T ) / variationCentralInstantPositionIoEarth[ 2 ];



//    /// w.r.t. position of second transmitter
//    variationInitialVectorEuropa = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionEuropaEarth[ 2 ] ).finished( );
//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) +  variationInitialVectorEuropa;
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );
//    modifiedrelativePositionEuropaSun = relativePositionEuropaSun - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaJupiter = relativePositionEuropaJupiter - variationInitialVectorEuropa;
//    modifiedrelativePositionEuropaIo = relativePositionEuropaIo - variationInitialVectorEuropa;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedrelativePositionJupiterEuropa = relativePositionJupiterEuropa + variationInitialVectorEuropa;
//    modifiedCentralInstantAccelerationJupiter = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( relativePositionJupiterIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationJupiter += calculatePointMassGravityAcceleration( modifiedrelativePositionJupiterEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= modifiedCentralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                                                                  rightAscensionAndDeclinationIo.first, rightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEuropaPosition[ 2 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 2 ];
//    numericalPartialsRwrtEuropaPosition[ 2 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 2 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEuropaPosition[ 2 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 2 ];


//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEuropaPosition[ 2 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 2 ];
//    numericalPartialsTwrtEuropaPosition[ 2 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 2 ];



//    /// w.r.t. position of receiver
//    variationInitialVectorEarth = ( Eigen::Vector3d( ) << 0.0, 0.0, variationCentralInstantPositionEuropaEarth[ 2 ] ).finished( );

//    modifiedInitialPositionEuropaEarth = ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedCentralInstantPositionIoEarth = ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ) -  variationInitialVectorEarth;
//    modifiedRightAscensionAndDeclinationIo = computeRightAscensionDeclination( modifiedCentralInstantPositionIoEarth );
//    modifiedrightAscensionAndDeclinationEuropa = computeRightAscensionDeclination( modifiedInitialPositionEuropaEarth );

//    modifiedRelativePositionEarthSun = relativePositionEarthSun - variationInitialVectorEarth;
//    modifiedRelativePositionEarthJupiter = relativePositionEarthJupiter - variationInitialVectorEarth;
//    modifiedRelativePositionEarthIo = relativePositionEarthIo - variationInitialVectorEarth;
//    modifiedRelativePositionEarthEuropa = relativePositionEarthEuropa - variationInitialVectorEarth;

//    modifiedRelativePositionIoEarth = relativePositionIoEarth + variationInitialVectorEarth;
//    modifiedrelativePositionEuropaEarth = relativePositionEuropaEarth + variationInitialVectorEarth;

//    modifiedcentralInstantAccelerationEarth = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEarth += calculatePointMassGravityAcceleration( modifiedRelativePositionEarthEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );

//    modifiedCentralInstantAccelerationIo = Eigen::Vector3d::Zero( );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( modifiedRelativePositionIoEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo += calculatePointMassGravityAcceleration( relativePositionIoEuropa, bodyMap[ "Europa" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedCentralInstantAccelerationIo -= centralInstantAccelerationJupiter;

//    modifiedcentralInstantAccelerationEuropa = Eigen::Vector3d::Zero( );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaSun, bodyMap[ "Sun" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( modifiedrelativePositionEuropaEarth, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaJupiter, bodyMap[ "Jupiter" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa += calculatePointMassGravityAcceleration( relativePositionEuropaIo, bodyMap[ "Io" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
//    modifiedcentralInstantAccelerationEuropa -= centralInstantAccelerationJupiter;

//    modifiedCentralInstantValueX = computeX( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );
//    modifiedcentralInstantValueY = computeY( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth );

//    modifiedFirstTimeDerivativeOfX = computePartialXwrtTime(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedFirstTimeDerivativeOfY = computePartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ) );

//    modifiedSecondTimeDerivativeOfX = computeSecondPartialXwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                                                                  modifiedRightAscensionAndDeclinationIo.first, modifiedRightAscensionAndDeclinationIo.second,
//                                                                  modifiedrightAscensionAndDeclinationEuropa.first, modifiedrightAscensionAndDeclinationEuropa.second );

//    modifiedSecondTimeDerivativeOfY = computeSecondPartialYwrtTime( modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                                                                  ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                                                                  modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth );

//    modifiedCubicPolynomialCoefficients = ( Eigen::Vector4d( ) << modifiedSecondTimeDerivativeOfX * modifiedSecondTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedSecondTimeDerivativeOfY,
//                                            3.0 * ( modifiedSecondTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedSecondTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedSecondTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedSecondTimeDerivativeOfY
//                                                    + modifiedFirstTimeDerivativeOfX * modifiedFirstTimeDerivativeOfX + modifiedFirstTimeDerivativeOfY * modifiedFirstTimeDerivativeOfY ),
//                                            2.0 * ( modifiedCentralInstantValueX * modifiedFirstTimeDerivativeOfX + modifiedcentralInstantValueY * modifiedFirstTimeDerivativeOfY ) ).finished( );

//    modifiedDepressedCubicPolynomialCoefficients = ( Eigen::Vector3d( ) << modifiedCubicPolynomialCoefficients[ 1 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 2 ] / modifiedCubicPolynomialCoefficients[ 0 ],
//            modifiedCubicPolynomialCoefficients[ 3 ] / modifiedCubicPolynomialCoefficients[ 0 ] ).finished( );

//    modifiedQ = computeConstantQ( modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );
//    modifiedR = computeConstantR( modifiedDepressedCubicPolynomialCoefficients[ 2 ], modifiedDepressedCubicPolynomialCoefficients[ 1 ], modifiedDepressedCubicPolynomialCoefficients[ 0 ] );

//    numericalPartialsQwrtEarthPosition[ 2 ] = ( modifiedQ - Q ) / variationCentralInstantPositionEuropaEarth[ 2 ];
//    numericalPartialsRwrtEarthPosition[ 2 ] = ( modifiedR - R ) / variationCentralInstantPositionEuropaEarth[ 2 ];

//    modifiedAngleTheta = computeAngleThetaSolutionCubicEquation(
//                modifiedCentralInstantPositionIoEarth, modifiedInitialPositionEuropaEarth,
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                modifiedCentralInstantAccelerationIo - modifiedcentralInstantAccelerationEarth, modifiedcentralInstantAccelerationEuropa - modifiedcentralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    numericalPartialsThetawrtEarthPosition[ 2 ] = ( modifiedAngleTheta - angleTheta ) / variationCentralInstantPositionEuropaEarth[ 2 ];


//    modifiedBeta = modifiedQ * modifiedQ * modifiedQ + modifiedR * modifiedR;
//    if ( ( modifiedR + sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedS = std::pow( modifiedR + sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedS = - std::pow( std::fabs( modifiedR + std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    if ( ( modifiedR - sqrt( modifiedBeta ) ) >= 0 )
//    {
//        modifiedT = std::pow( modifiedR - sqrt( modifiedBeta ), 1.0 / 3.0 );
//    }
//    else
//    {
//        modifiedT = - std::pow( std::fabs( modifiedR - std::sqrt( modifiedBeta ) ), 1.0 / 3.0 );
//    }
//    numericalPartialsSwrtEarthPosition[ 2 ] = ( modifiedS - S ) / variationCentralInstantPositionEuropaEarth[ 2 ];
//    numericalPartialsTwrtEarthPosition[ 2 ] = ( modifiedT - T ) / variationCentralInstantPositionEuropaEarth[ 2 ];


//    //////////////////////////////////////////////////////////////////////////////////



//    Eigen::Vector3d partialsQWrtIoPosition = computePartialsQWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    std::cout << "partials intermediate Q w.r.t. first transmitter position (function): " << partialsQWrtIoPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate Q w.r.t. first transmitter position: " << numericalPartialsQwrtIoPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsRWrtIoPosition = computePartialsRWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    std::cout << "partials intermediate R w.r.t. first transmitter position (function): " << partialsRWrtIoPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate R w.r.t. first transmitter position: " << numericalPartialsRwrtIoPosition.transpose( ) << "\n\n";



//    Eigen::Vector3d partialsQWrtEuropaPosition = computePartialsQWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false );

//    std::cout << "partials intermediate Q w.r.t. second transmitter position (function): " << partialsQWrtEuropaPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate Q w.r.t. second transmitter position: " << numericalPartialsQwrtEuropaPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsRWrtEuropaPosition = computePartialsRWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false );

//    std::cout << "partials intermediate R w.r.t. second transmitter position (function): " << partialsRWrtEuropaPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate R w.r.t. second transmitter position: " << numericalPartialsRwrtEuropaPosition.transpose( ) << "\n\n";


//    Eigen::Vector3d partialsQWrtEarthPosition = - computePartialsQWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true )
//            - computePartialsQWrtCartesianPosition(
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                            centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                            - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );

//    std::cout << "partials intermediate Q w.r.t. receiver position (function): " << partialsQWrtEarthPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate Q w.r.t. receiver position: " << numericalPartialsQwrtEarthPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsRWrtEarthPosition = - computePartialsRWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true )
//            - computePartialsRWrtCartesianPosition(
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                            centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                            - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );

//    std::cout << "partials intermediate R w.r.t. receiver position (function): " << partialsRWrtEarthPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate R w.r.t. receiver position: " << numericalPartialsRwrtEarthPosition.transpose( ) << "\n\n";


////    Eigen::Vector3d partialsAngleThetaWrtIoPosition = computePartialsThetaWrtCartesianPosition(
////                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
////                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
////                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
////                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );
////    std::cout << "partials angle theta cubic equation w.r.t. first transmitter position (function): " << partialsAngleThetaWrtIoPosition.transpose( ) << "\n\n";
////    std::cout << "numerical partials angle theta cubic equation w.r.t. first transmitter position: " << numericalPartialsThetawrtIoPosition.transpose( ) << "\n\n";

////    Eigen::Vector3d partialsAngleThetaWrtEuropaPosition = computePartialsThetaWrtCartesianPosition(
////                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
////                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
////                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
////                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false );
////    std::cout << "partials angle theta cubic equation w.r.t. second transmitter position (function): " << partialsAngleThetaWrtEuropaPosition.transpose( ) << "\n\n";
////    std::cout << "numerical partials angle theta cubic equation w.r.t. second transmitter position: " << numericalPartialsThetawrtEuropaPosition.transpose( ) << "\n\n";

////    Eigen::Vector3d partialsAngleThetaWrtEarthPosition = - computePartialsThetaWrtCartesianPosition(
////                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
////                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
////                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
////                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true )
////            - computePartialsThetaWrtCartesianPosition(
////                            ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
////                            ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
////                            centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
////                            - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );
////    std::cout << "partials angle theta cubic equation w.r.t. receiver position (function): " << partialsAngleThetaWrtEarthPosition.transpose( ) << "\n\n";
////    std::cout << "numerical partials angle theta cubic equation w.r.t. receiver position: " << numericalPartialsThetawrtEarthPosition.transpose( ) << "\n\n";


//    Eigen::Vector3d partialsTWrtIoPosition = computePartialsTWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true, T );
//    std::cout << "partials intermediate variable T w.r.t. first transmitter position (function): " << partialsTWrtIoPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate T w.r.t. first transmitter position: " << numericalPartialsTwrtIoPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsSWrtIoPosition = computePartialsSWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true, S );
//    std::cout << "partials intermediate variable S w.r.t. first transmitter position (function): " << partialsSWrtIoPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate S w.r.t. first transmitter position: " << numericalPartialsSwrtIoPosition.transpose( ) << "\n\n";


//    Eigen::Vector3d partialsTWrtEuropaPosition = computePartialsTWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false, T );
//    std::cout << "partials intermediate variable T w.r.t. second transmitter position (function): " << partialsTWrtEuropaPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate T w.r.t. second transmitter position: " << numericalPartialsTwrtEuropaPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsSWrtEuropaPosition = computePartialsSWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false, S );
//    std::cout << "partials intermediate variable S w.r.t. second transmitter position (function): " << partialsSWrtEuropaPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate S w.r.t. second transmitter position: " << numericalPartialsSwrtEuropaPosition.transpose( ) << "\n\n";


//    Eigen::Vector3d partialsTWrtEarthPosition = - computePartialsTWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true, T )
//            - computePartialsTWrtCartesianPosition(
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                            centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                            - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false, T );
//    std::cout << "partials intermediate variable T w.r.t. receiver position (function): " << partialsTWrtEarthPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate T w.r.t. receiver position: " << numericalPartialsTwrtEarthPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsSWrtEarthPosition = - computePartialsSWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true, S )
//            - computePartialsSWrtCartesianPosition(
//                           ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                           ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                           centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                           - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false, S );
//    std::cout << "partials intermediate variable S w.r.t. receiver position (function): " << partialsSWrtEarthPosition.transpose( ) << "\n\n";
//    std::cout << "numerical partials intermediate S w.r.t. receiver position: " << numericalPartialsSwrtEarthPosition.transpose( ) << "\n\n";




//    Eigen::Vector3d partialsCentralInstantWrtIoPosition = computePartialsCentralInstantWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, true );

//    std::cout << "partials central instant w.r.t. first transmitter position (function): " << partialsCentralInstantWrtIoPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsCentralInstantWrtEuropaPosition = computePartialsCentralInstantWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                partialsAccelerationIoWrtPositionIo, partialscentralInstantAccelerationEuropa, false );

//    std::cout << "partials central instant w.r.t. second transmitter position (function): " << partialsCentralInstantWrtEuropaPosition.transpose( ) << "\n\n";

//    Eigen::Vector3d partialsCentralInstantWrtEarthPosition = - computePartialsCentralInstantWrtCartesianPosition(
//                ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, true )
//            - computePartialsCentralInstantWrtCartesianPosition(
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 0, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 0, 3 ),
//                            ( centralInstantIoState - centralInstantEarthState ).segment( 3, 3 ), ( centralInstantEuropaState - centralInstantEarthState ).segment( 3, 3 ),
//                            centralInstantAccelerationIo - centralInstantAccelerationEarth, centralInstantAccelerationEuropa - centralInstantAccelerationEarth,
//                            - partialsAccelerationIoWrtPositionEarth, - partialsAccelerationEuropaWrtPositionEarth, false );

//    std::cout << "partials central instant w.r.t. receiver position (function): " << partialsCentralInstantWrtEarthPosition.transpose( ) << "\n\n";


//    std::cout << "initial value central instant: " << estimatedCentralInstant << "\n\n";
//    std::cout << "estimated central instant: " << - 1.0 / 3.0 * depressedCubicPolynomialCoefficients[ 0 ] + ( S + T ) << "\n\n";


////    std::cout << "TEST TEST TEST: " << centralInstantValueX * firstTimePartialOfX + centralInstantValueY * firstTimePartialOfY << "\n\n";

////    double newEstimatedCentralInstant = computeCentralInstantFromDependentVariables( simulationStartEpoch, 60.0,
////                                                                                     modifiedDependentVariablesHistory );
////    std::cout << "newEstimatedCentralInstant: " << newEstimatedCentralInstant << "\n\n";


    /////////////////////////////////////////////////////////////////////////////////////////


//    std::cout.precision( 20 );

////    // Define and create ground stations.
////    std::vector< std::pair< std::string, std::string > > groundStations;
////    groundStations.resize( 2 );
////    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
////    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


//    // Initial guess central instant.
//    double estimatedCentralInstant = 116200.0;

//    // Specify initial time
//    double initialEphemerisTime = 0.0; //estimatedCentralInstant - 3600.0;
//    double finalEphemerisTime = 1.0 * physical_constants::JULIAN_YEAR; // estimatedCentralInstant + 3600.0;

//    // Load spice kernel.
//    spice_interface::loadStandardSpiceKernels( );

//    // Define body settings for simulation.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );
//    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Jupiter" );
//    bodiesToCreate.push_back( "Io" );
//    bodiesToCreate.push_back( "Europa" );

//    // Create body objects.
//    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
//            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime, finalEphemerisTime );
//    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
//    {
//        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
//        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
//    }
//    NamedBodyMap bodyMap = createBodies( bodySettings );

////    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
////    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
////                "Earth", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
////    bodyMap[ "Earth" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );
////    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
////                "Io", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
////    bodyMap[ "Io" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );
////    bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
////                "Europa", "SSB", "ECLIPJ2000", "NONE", 1.1e7 );
////    bodyMap[ "Europa" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "J2000" ) );

//    // Finalize body creation.
//    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

//    SelectedAccelerationMap singleArcAccelerationMap;
//    std::vector< std::string > singleArcBodiesToPropagate = { "Earth", "Io", "Europa" };
//    std::vector< std::string > singleArcCentralBodies = { "Sun", "Jupiter", "Jupiter" };

//    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
//    {
//        singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ "Jupiter" ].push_back(
//                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
//        singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ "Sun" ].push_back(
//                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

//        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
//        {
//            if( i != j )
//            {
//                singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ singleArcBodiesToPropagate.at( j ) ].push_back(
//                            std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
//            }
//        }
//    }

//    basic_astrodynamics::AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
//                bodyMap, singleArcAccelerationMap, singleArcBodiesToPropagate, singleArcCentralBodies );

//    Eigen::VectorXd singleArcInitialState = propagators::getInitialStatesOfBodies(
//                singleArcBodiesToPropagate, singleArcCentralBodies, bodyMap, initialEphemerisTime );

//    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >
//                                      ( propagators::total_acceleration_dependent_variable, "Io" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Io", "Io", "Jupiter" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Io", "Earth", "Jupiter" ) );
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >
//                                      ( propagators::total_acceleration_dependent_variable, "Europa" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Europa", "Europa", "Jupiter" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Europa", "Earth", "Jupiter" ) );
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >
//                                      ( propagators::total_acceleration_dependent_variable, "Earth" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Earth", "Earth", "Sun" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Earth", "Io", "Sun" ) );
//    dependentVariablesList.push_back(
//                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
//                    "Earth", "Europa", "Sun" ) );

//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >
//                                      ( propagators::relative_position_dependent_variable, "Io", "Earth" ) );
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >
//                                      ( propagators::relative_position_dependent_variable, "Europa", "Earth" ) );

//    dependentVariablesList.push_back(
//                std::make_shared< propagators::AccelerationPartialWrtStateSaveSettings >(
//                    "Io", "Jupiter", basic_astrodynamics::point_mass_gravity, "Io", "Jupiter" ) );

//    // Create object with list of dependent variables
//    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
//            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings =
//            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//            ( singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToPropagate, singleArcInitialState, finalEphemerisTime,
//              propagators::cowell, dependentVariablesToSave );

//    const double fixedStepSize = 3600.0;
//    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
//            std::make_shared< numerical_integrators::IntegratorSettings< > >( numerical_integrators::rungeKutta4, 0.0, fixedStepSize );

//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
//    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatableVectorParameters;
//    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedInitialStateParameters;

//    estimatedInitialStateParameters.push_back(
//                std::make_shared< InitialTranslationalStateParameter< double > >(
//                    "Earth", propagators::getInitialStateOfBody( "Earth", "SSB", bodyMap, initialEphemerisTime ) ) );
//    estimatedInitialStateParameters.push_back(
//                std::make_shared< InitialTranslationalStateParameter< double > >(
//                    "Io", propagators::getInitialStateOfBody( "Io", "SSB", bodyMap, initialEphemerisTime ) ) );
//    estimatedInitialStateParameters.push_back(
//                std::make_shared< InitialTranslationalStateParameter< double > >(
//                    "Europa", propagators::getInitialStateOfBody( "Europa", "SSB", bodyMap, initialEphemerisTime ) ) );

//    std::shared_ptr< EstimatableParameter< double > > relativisticParameter;
//    relativisticParameter = std::make_shared< estimatable_parameters::PPNParameterGamma >( );
////        relativisticParameter = std::make_shared< estimatable_parameters::EquivalencePrincipleLpiViolationParameter >( );

//    estimatableDoubleParameters.push_back( relativisticParameter );

//    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
//            std::make_shared< EstimatableParameterSet< double > >( estimatableDoubleParameters,
//                    estimatableVectorParameters, estimatedInitialStateParameters );

////    // Define list of parameters to estimate.
////    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
////    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
////                                  "AlienSpaceship", alienSpaceshipInitialState, "Phobos" ) );
////    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "AlienSpaceship", radiation_pressure_coefficient ) );
////    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
////    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
////                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
////    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
////                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

////    // Create parameters
////    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
////            createParametersToEstimate( parameterNames, bodyMap );

////    // Print identifiers and indices of parameters to terminal.
////    printEstimatableParameterEntries( parametersToEstimate );

//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create simulation object and propagate dynamics.
//    propagators::SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
//                bodyMap, integratorSettings, singleArcPropagatorSettings, parametersToEstimate, true,
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true,
//                dependentVariablesToSave );

//    std::map< double, Eigen::VectorXd > dependentVariablesHistory = variationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );
//    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ; itr!= dependentVariablesHistory.end( ) ; itr++ )
//    {
//        Eigen::Matrix< double, 3, 1 > sphericalCoordinatesIoSeenFromEarth = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                    itr->second.segment( 63, 3 ) ).template cast< double >( );
//        double rightAscensionIo = sphericalCoordinatesIoSeenFromEarth.z( );
//        double declinationIo = mathematical_constants::PI / 2.0 - sphericalCoordinatesIoSeenFromEarth.y( );

//        Eigen::Matrix< double, 3, 1 > sphericalCoordinatesEuropaSeenFromEarth = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
//                    itr->second.segment( 66, 3 ) ).template cast< double >( );
//        double rightAscensionEuropa = sphericalCoordinatesEuropaSeenFromEarth.z( );
//        double declinationEuropa = mathematical_constants::PI / 2.0 - sphericalCoordinatesEuropaSeenFromEarth.y( );

//        double averageDeclination = ( declinationIo + declinationEuropa ) / 2.0;

//        double apparentDistance = std::sqrt( ( rightAscensionIo - rightAscensionEuropa ) * cos( averageDeclination )
//                                             * ( rightAscensionIo - rightAscensionEuropa ) * cos( averageDeclination )
//                                             + ( declinationIo - declinationEuropa ) * ( declinationIo - declinationEuropa ) );

//        std::cout << "time: " << itr->first << " ; apparent distance Io-Europa: " << apparentDistance * 180.0 / mathematical_constants::PI * 3600.0 << " arcseconds" << "\n\n";
//    }

//    // Retrieve dependent variables interface.
//    std::shared_ptr< propagators::SingleArcDependentVariablesInterface > dependentVariablesInterface =
//            std::dynamic_pointer_cast< propagators::SingleArcDependentVariablesInterface >(
//            variationalEquationsSimulator.getDependentVariablesInterface( ) );

//    // Total acceleration dependent variable settings.
//    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
//            = std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                propagators::total_acceleration_dependent_variable, "Io" );

//    // Partial of total acceleration w.r.t. translational state dependent variable settings.
//    std::shared_ptr< propagators::TotalAccelerationPartialWrtStateSaveSettings > partialTotalAccelerationWrtStateDependentVariable
//            = std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >( "Io", "Io", "Jupiter" );

//    std::cout << "total acceleration from interface: " << dependentVariablesInterface->getSingleDependentVariable(
//                     totalAccelerationDependentVariable, ( initialEphemerisTime - finalEphemerisTime ) / 2.0 ).transpose( ) << "\n\n";

//    std::cout << "partial of total acceleration from interface: " << dependentVariablesInterface->getSingleDependentVariable(
//                     partialTotalAccelerationWrtStateDependentVariable, ( initialEphemerisTime - finalEphemerisTime ) / 2.0 ).transpose( ) << "\n\n";

//    // Io-Earth relative position dependent variable settings.
//    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > ioEarthRelativePositionDependentVariable
//            = std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                propagators::relative_position_dependent_variable, "Io", "Earth" );

//    // Europa-Earth relative position dependent variable settings.
//    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > europaEarthRelativePositionDependentVariable
//            = std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                propagators::relative_position_dependent_variable, "Europa", "Earth" );

//    std::cout << "Io-Earth relative position from interface: " << dependentVariablesInterface->getSingleDependentVariable(
//                     ioEarthRelativePositionDependentVariable, 31168800.0 /*( initialEphemerisTime - finalEphemerisTime ) / 2.0*/ ).transpose( ) << "\n\n";

//    std::cout << "Europa-Earth relative position from interface: " << dependentVariablesInterface->getSingleDependentVariable(
//                     europaEarthRelativePositionDependentVariable, 31168800.0 /*( initialEphemerisTime - finalEphemerisTime ) / 2.0*/ ).transpose( ) << "\n\n";






////    basic_astrodynamics::AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
////                bodyMap, singleArcAccelerationMap, singleArcBodiesToPropagate, singleArcCentralBodies );

////    Eigen::VectorXd singleArcInitialState = getInitialStatesOfBodies(
////                singleArcBodiesToPropagate, singleArcCentralBodies, bodyMap, initialTime );

//    // Define link ends for observations.
//    LinkEnds linkEnds;
//    linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );
//    linkEnds[ transmitter ] = std::make_pair( "Io" , ""  );
//    linkEnds[ transmitter2 ] = std::make_pair( "Europa" , ""  );



//    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
////    parameterPerturbationMultipliers( 2 ) = 10.0;
//    // Test partials with constant ephemerides (allows test of position partials)
//    {
////        // Create environment
////        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

////        // Set link ends for observation model
////        LinkEnds linkEnds;
////        linkEnds[ transmitter ] = groundStations[ 1 ];
////        linkEnds[ receiver ] = groundStations[ 0 ];


//        // Create light-time correction settings
//        std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
//        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
//        lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
//                                                    lightTimePerturbingBodies ) );

//        // Create observation settings
//        std::shared_ptr< ObservationSettings > observableSettings = std::make_shared< ObservationSettings >
//                ( mutual_approximation, lightTimeCorrectionSettings, std::make_shared< ConstantObservationBiasSettings >(
//                      ( Eigen::Vector1d( ) << 1.0 / 3600.0 * mathematical_constants::PI / 180.0 ).finished( ), true ) );

//        // Create mutual approximation model.
//        std::shared_ptr< ObservationModel< 1, double, double > > mutualApproximationModel =
//               ObservationModelCreator< 1, double, double >::createObservationModel(
//                    linkEnds, observableSettings, bodyMap );
////        std::shared_ptr< ObservationBias< 1 > > observationBias = mutualApproximationModel->getObservationBiasCalculator( );

////        std::vector< double > linkEndTimes;
////        std::vector< Eigen::Vector6d > linkEndStates;
////        Eigen::Vector1d observationFromReceptionTime = mutualApproximationModel->computeObservations( estimatedCentralInstant, receiver );

////        // Generate mutual approximation model
////        std::vector< std::string > perturbingBodies;
////        perturbingBodies.push_back( "Sun" );
////        std::shared_ptr< ObservationModel< 1 > > mutualApproximationModel =
////                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
////                    linkEnds, std::make_shared< observation_models::ObservationSettings >(
////                        observation_models::mutual_approximation, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
////                            perturbingBodies ) ), bodyMap  );

////        std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
////        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatableVectorParameters;
////        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedInitialStateParameters;

////        estimatedInitialStateParameters.push_back(
////                    std::make_shared< InitialTranslationalStateParameter< double > >(
////                        "Earth", propagators::getInitialStateOfBody(
////                            "Earth", "SSB", bodyMap, 1.1e7 ) ) );
////        estimatedInitialStateParameters.push_back(
////                    std::make_shared< InitialTranslationalStateParameter< double > >(
////                        "Io", propagators::getInitialStateOfBody(
////                            "Io", "SSB", bodyMap, 1.1e7 ) ) );
////        estimatedInitialStateParameters.push_back(
////                    std::make_shared< InitialTranslationalStateParameter< double > >(
////                        "Europa", propagators::getInitialStateOfBody(
////                            "Europa", "SSB", bodyMap, 1.1e7 ) ) );

////        std::shared_ptr< EstimatableParameter< double > > relativisticParameter;
////        relativisticParameter = std::make_shared< estimatable_parameters::PPNParameterGamma >( );
//////        relativisticParameter = std::make_shared< estimatable_parameters::EquivalencePrincipleLpiViolationParameter >( );

////        estimatableDoubleParameters.push_back( relativisticParameter );

////        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
////                std::make_shared< EstimatableParameterSet< double > >( estimatableDoubleParameters,
////                        estimatableVectorParameters, estimatedInitialStateParameters );



//////        // Create parameter objects.
//////        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
//////                createEstimatableParameters( bodyMap, 1.1E7 );


//        /////////// STUFF TAKEN FROM testObservationPartials FUNCTION

//        // Create observation partials.
//        std::map< LinkEnds, std::shared_ptr< ObservationModel< 1 > > > observationModelList;
//        observationModelList[ linkEnds ] = mutualApproximationModel;

//        std::shared_ptr< ObservationPartialCreator< 1, double, double > > observationPartialCreator =
//                std::make_shared< ObservationPartialCreator< 1, double, double > >( );
//        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
//                std::shared_ptr< PositionPartialScaling > > fullAnalyticalPartialSet =
//                observationPartialCreator->createObservationPartials(
//                    mutual_approximation, observationModelList, bodyMap, parametersToEstimate, dependentVariablesInterface ).begin( )->second;
//        std::shared_ptr< PositionPartialScaling > positionPartialScaler = fullAnalyticalPartialSet.second;

//        // Iterate over link ends, compute and test partials for observable referenced at each link end.
//        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
//             linkEndIterator++ )
//        {

//            if ( linkEndIterator->first == receiver )
//            {

//                std::cout << "============================ REFERENCE LINK END =============" << linkEndIterator->first << std::endl;
//                // Evaluate nominal observation values
//                std::vector< Eigen::Vector6d > vectorOfStates;
//                std::vector< double > vectorOfTimes;
//                double observationTime = 31168800.0;

//                Eigen::VectorXd currentObservation = mutualApproximationModel->computeObservationsWithLinkEndData(
//                            observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates );

//                // Calculate analytical observation partials.
//                if( positionPartialScaler != NULL )
//                {
//                    positionPartialScaler->update( vectorOfStates, vectorOfTimes, static_cast< LinkEndType >( linkEndIterator->first ), linkEnds,
//                                                   currentObservation );
//                }

//                typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > >
//                        ObservationPartialReturnType;
//                std::vector< ObservationPartialReturnType > analyticalObservationPartials =
//                        calculateAnalyticalPartials< 1 >(
//                            fullAnalyticalPartialSet.first, vectorOfStates, vectorOfTimes, linkEndIterator->first, currentObservation );

//            }

//        }

//        std::cout << "dependent variables from interface at central instant t0: " <<
//                     dependentVariablesInterface->getDependentVariables( 31168800.0 ).transpose( ) << "\n\n";

//        std::cout << "partial of total acceleration of Io w.r.t. its own state from dependent variables (using setting): " <<
//                     dependentVariablesInterface->getSingleDependentVariable(
//                         std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >( "Io", "Io", "Jupiter" ), 31168800.0 ).transpose( ) << "\n\n";

//        std::cout << "partial of total acceleration of Io w.r.t. its own state from dependent variables (using ID): " <<
//                     dependentVariablesInterface->getSingleDependentVariable(
//                         propagators::getDependentVariableId( std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >( "Io", "Io", "Jupiter" ) ),
//                         18, 31168800.0 ).transpose( ) << "\n\n";


//        /////////////////////////// END OF STUFF TAKEN FROM testObservationPartials FUNCTION



////        testObservationPartials( mutualApproximationModel, bodyMap, parametersToEstimate, linkEnds,
////                                 mutual_approximation, 1.0E-4, true, true, 1.0, parameterPerturbationMultipliers, dependentVariablesInterface, 31168800.0 );
//    }


////    // Test partials with real ephemerides (without test of position partials)
////    {
////        std::cout << "Test 1" << std::endl;
////        // Create environment
////        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

////        // Set link ends for observation model
////        LinkEnds linkEnds;
////        linkEnds[ transmitter ] = groundStations[ 1 ];
////        linkEnds[ receiver ] = groundStations[ 0 ];

////        // Generate one-way range model
////        std::shared_ptr< ObservationModel< 2 > > angularPositionModel =
////                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
////                    linkEnds, std::make_shared< observation_models::ObservationSettings >(
////                        observation_models::angular_position ), bodyMap  );

////        // Create parameter objects.
////        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
////                createEstimatableParameters( bodyMap, 1.1E7 );

////        testObservationPartials( angularPositionModel, bodyMap, fullEstimatableParameterSet, linkEnds, angular_position, 1.0E-4,
////        false, true, 1.0, parameterPerturbationMultipliers );

////    }
}



BOOST_AUTO_TEST_CASE( testCentralInstantPartials )
{
    std::cout.precision( 16 );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 1300000.0; // 41680.0; //0.0;
    const double simulationEndEpoch = 2.0 * tudat::physical_constants::JULIAN_YEAR;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Io" );
    bodiesToCreate.push_back( "Europa" );
    bodiesToCreate.push_back( "Ganymede" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsEarth;
    currentAccelerationsEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEarth[ "Ganymede" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsJupiter;
    currentAccelerationsJupiter[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsJupiter[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsJupiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsJupiter[ "Ganymede" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsIo;
    currentAccelerationsIo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsIo[ "Ganymede" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsEuropa;
    currentAccelerationsEuropa[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsEuropa[ "Ganymede" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerationsGanymede;
    currentAccelerationsGanymede[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsGanymede[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsGanymede[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsGanymede[ "Io" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    currentAccelerationsGanymede[ "Europa" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );


    accelerationMap[ "Ganymede" ] = currentAccelerationsGanymede;
    accelerationMap[ "Jupiter" ] = currentAccelerationsJupiter;
    accelerationMap[ "Io" ] = currentAccelerationsIo;
    accelerationMap[ "Europa" ] = currentAccelerationsEuropa;
    accelerationMap[ "Earth" ] = currentAccelerationsEarth;

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Earth" );
    bodiesToPropagate.push_back( "Jupiter" );
    bodiesToPropagate.push_back( "Io" );
    bodiesToPropagate.push_back( "Europa" );
    bodiesToPropagate.push_back( "Ganymede" );

    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "SSB" );
    centralBodies.push_back( "SSB" );
    centralBodies.push_back( "Jupiter" );
    centralBodies.push_back( "Jupiter" );
    centralBodies.push_back( "Jupiter" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, simulationStartEpoch );
//    systemInitialState[ 7 ] *= 0.999;


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

    double limitValueImpactParameter = 10.0 / 3600.0 * mathematical_constants::PI / 180.0;
    std::cout << "limit value impact parameter: " << limitValueImpactParameter << "\n\n";


    for ( unsigned int currentBody = 2 ; currentBody < bodiesToPropagate.size( ) ; currentBody++ )
    {
        dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                        relative_position_dependent_variable, bodiesToPropagate[ currentBody ], "Earth" ) );
    }

    for ( unsigned int currentBody = 2 ; currentBody < bodiesToPropagate.size( ) - 1 ; currentBody++ )
    {

        for ( unsigned int otherBody = currentBody + 1 ; otherBody < bodiesToPropagate.size( ) ; otherBody ++ )
        {
            dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                  relative_distance_dependent_variable, bodiesToPropagate[ currentBody ],
                                                  bodiesToPropagate[ otherBody ] ) );
        }
    }

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > secondDependentVariablesList = dependentVariablesList;
    secondDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_velocity_dependent_variable, "Io", "Earth" ) );
    secondDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_velocity_dependent_variable, "Europa", "Earth" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( secondDependentVariablesList );


    std::string observer = "Earth";
    std::vector< std::string > listOfObjects;
    listOfObjects.push_back( "Io" );
    listOfObjects.push_back( "Europa" );
    listOfObjects.push_back( "Ganymede" );
    bool useThresholdAsLowerLimit = true;;
//    std::pair< std::string, std::string > bodiesInvolvedInMutualApproximations;
    std::function< bool( const double ) > customTerminationFunction =
            std::bind( &isApparentDistanceBelowThreshold, std::placeholders::_1, limitValueImpactParameter,
                       "Earth", listOfObjects, bodyMap, simulationEndEpoch, useThresholdAsLowerLimit); //, bodiesInvolvedInMutualApproximations );

    std::shared_ptr< PropagationTerminationSettings > customTerminationSettings =
            std::make_shared< PropagationCustomTerminationSettings >( customTerminationFunction );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, customTerminationSettings, cowell, dependentVariablesToSave );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, 3600.0 / 10.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBITS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory( );

    double epochEndFirstPropagation = integrationResult.rbegin( )->first;

    Eigen::VectorXd finalPropagatedState = integrationResult.rbegin( )->second;
    Eigen::VectorXd finalPropagatedStateIo = finalPropagatedState.segment( 12, 6 );
    Eigen::VectorXd finalPropagatedStateEuropa = finalPropagatedState.segment( 18, 6 );
    Eigen::VectorXd finalPropagatedStateGanymede = finalPropagatedState.segment( 24, 6 );

    Eigen::VectorXd finalDependentVariablesValues = dependentVariablesHistory.rbegin( )->second;

    double finalDistanceIoEuropa = finalDependentVariablesValues[ 9 ];
    double finalDistanceIoGanymede = finalDependentVariablesValues[ 10 ];
    double finalDistanceEuropaGanymede = finalDependentVariablesValues[ 11 ];

    Eigen::VectorXd finalRelativePositionIoEarth = finalDependentVariablesValues.segment( 0, 3 );
    Eigen::Matrix< double, 3, 1 > finalSphericalStateIo = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                finalRelativePositionIoEarth.segment( 0, 3 ) ).template cast< double >( );
    double rightAscensionIo = finalSphericalStateIo.z( );
    double declinationIo = mathematical_constants::PI / 2.0 - finalSphericalStateIo.y( );

    Eigen::VectorXd finalRelativePositionEuropaEarth = finalDependentVariablesValues.segment( 3, 3 );
    Eigen::Matrix< double, 3, 1 > finalSphericalStateEuropa = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                finalRelativePositionEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
    double rightAscensionEuropa = finalSphericalStateEuropa.z( );
    double declinationEuropa = mathematical_constants::PI / 2.0 - finalSphericalStateEuropa.y( );

    Eigen::VectorXd finalRelativePositionGanymedeEarth = finalDependentVariablesValues.segment( 6, 3 );
    Eigen::Matrix< double, 3, 1 > finalSphericalStateGanymede = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                finalRelativePositionGanymedeEarth.segment( 0, 3 ) ).template cast< double >( );
    double rightAscensionGanymede = finalSphericalStateGanymede.z( );
    double declinationGanymede = mathematical_constants::PI / 2.0 - finalSphericalStateGanymede.y( );

    double deltaRightAscensionIoEuropa = rightAscensionEuropa - rightAscensionIo;
    double deltaDeclinationIoEuropa = declinationEuropa - declinationIo;
    double apparentDistanceIoEuropa = std::sqrt( ( deltaRightAscensionIoEuropa * std::cos( declinationIo ) ) * ( deltaRightAscensionIoEuropa * std::cos( declinationIo ) )
                                                 + deltaDeclinationIoEuropa * deltaDeclinationIoEuropa );

    double deltaRightAscensionIoGanymede = rightAscensionGanymede - rightAscensionIo;
    double deltaDeclinationIoGanymede = declinationGanymede - declinationIo;
    double apparentDistanceIoGanymede = std::sqrt( ( deltaRightAscensionIoGanymede * std::cos( declinationIo ) ) * ( deltaRightAscensionIoGanymede * std::cos( declinationIo ) )
                                                 + deltaDeclinationIoGanymede * deltaDeclinationIoGanymede );

    double deltaRightAscensionEuropaGanymede = rightAscensionGanymede - rightAscensionEuropa;
    double deltaDeclinationEuropaGanymede = declinationGanymede - declinationEuropa;
    double apparentDistanceEuropaGanymede = std::sqrt( ( deltaRightAscensionEuropaGanymede * std::cos( declinationEuropa ) ) * ( deltaRightAscensionEuropaGanymede * std::cos( declinationEuropa ) )
                                                 + deltaDeclinationEuropaGanymede * deltaDeclinationEuropaGanymede );

    std::cout << "apparent distance Io-Europa: " << apparentDistanceIoEuropa * 180.0 / mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
    std::cout << "apparent distance Io-Ganymede: " << apparentDistanceIoGanymede * 180.0 / mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";
    std::cout << "apparent distance Europa-Ganymede: " << apparentDistanceEuropaGanymede * 180.0 / mathematical_constants::PI * 3600.0 << " arcseconds." << "\n\n";


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////           SECOND PROPAGATION CLOSE MUTUAL APPROXIMATION          ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd newSystemInitialState = finalPropagatedState;


    useThresholdAsLowerLimit = false;
    std::vector< std::string > newListOfObjects;
    newListOfObjects.push_back( "Io" );
    newListOfObjects.push_back( "Europa" );
    std::function< bool( const double ) > newCustomTerminationFunction =
            std::bind( &isApparentDistanceBelowThreshold, std::placeholders::_1, limitValueImpactParameter,
                       "Earth", newListOfObjects, bodyMap, simulationEndEpoch, useThresholdAsLowerLimit );

    std::shared_ptr< PropagationTerminationSettings > newCustomTerminationSettings =
            std::make_shared< PropagationCustomTerminationSettings >( newCustomTerminationFunction );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > newPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, newSystemInitialState,
                newCustomTerminationSettings, cowell, dependentVariablesToSave );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > newIntegratorSettings = std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, epochEndFirstPropagation, 1.0 );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > newDynamicsSimulator( bodyMap, newIntegratorSettings, newPropagatorSettings, true, false, true );

    std::map< double, Eigen::VectorXd > newIntegrationResult = newDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > newDependentVariablesHistory = newDynamicsSimulator.getDependentVariableHistory( );


    // Create map with apparent relative distances.
    std::map< double, Eigen::VectorXd > apparentDistancesHistory;
    std::map< double, Eigen::VectorXd > angularPositionsHistory;

    std::map< double, Eigen::VectorXd > XYhistory;

    for ( std::map< double, Eigen::VectorXd >::iterator itr = newDependentVariablesHistory.begin( ) ;
          itr != newDependentVariablesHistory.end( ) ; itr++ )
    {
        Eigen::VectorXd relativePositionIoEarth = itr->second.segment( 0, 3 );
        Eigen::VectorXd relativePositionEuropaEarth = itr->second.segment( 3, 3 );
        Eigen::VectorXd relativePositionGanymedeEarth = itr->second.segment( 6, 3 );

        Eigen::VectorXd relativeVelocityIoEarth = itr->second.segment( 12, 3 );
        Eigen::VectorXd relativeVelocityEuropaEarth = itr->second.segment( 15, 3 );

        Eigen::Matrix< double, 3, 1 > sphericalStateIo = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                    relativePositionIoEarth.segment( 0, 3 ) ).template cast< double >( );
        double currentRightAscensionIo = sphericalStateIo.z( );
        double currentDeclinationIo = mathematical_constants::PI / 2.0 - sphericalStateIo.y( );

        Eigen::Matrix< double, 3, 1 > sphericalStateEuropa = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                    relativePositionEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
        double currentRightAscensionEuropa = sphericalStateEuropa.z( );
        double currentDeclinationEuropa = mathematical_constants::PI / 2.0 - sphericalStateEuropa.y( );

        Eigen::Matrix< double, 3, 1 > sphericalStateGanymede = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                    relativePositionGanymedeEarth.segment( 0, 3 ) ).template cast< double >( );
        double currentRightAscensionGanymede = sphericalStateGanymede.z( );
        double currentDeclinationGanymede = mathematical_constants::PI / 2.0 - sphericalStateGanymede.y( );

        double differenceRightAscensionIoEuropa = currentRightAscensionEuropa - currentRightAscensionIo;
        double differenceDeclinationIoEuropa = currentDeclinationEuropa - currentDeclinationIo;
        double currentApparentDistanceIoEuropa = std::sqrt( ( differenceRightAscensionIoEuropa * std::cos( currentDeclinationIo ) )
                                                            * ( differenceRightAscensionIoEuropa * std::cos( currentDeclinationIo ) )
                                                     + differenceDeclinationIoEuropa * differenceDeclinationIoEuropa );

        double differenceRightAscensionIoGanymede = currentRightAscensionGanymede - currentRightAscensionIo;
        double differenceDeclinationIoGanymede = currentDeclinationGanymede - currentDeclinationIo;
        double currentApparentDistanceIoGanymede = std::sqrt( ( differenceRightAscensionIoGanymede * std::cos( currentDeclinationIo ) )
                                                              * ( differenceRightAscensionIoGanymede * std::cos( currentDeclinationIo ) )
                                                     + differenceDeclinationIoGanymede * differenceDeclinationIoGanymede );

        double differenceRightAscensionEuropaGanymede = currentRightAscensionGanymede - currentRightAscensionEuropa;
        double differenceDeclinationEuropaGanymede = currentDeclinationGanymede - currentDeclinationEuropa;
        double currentApparentDistanceEuropaGanymede = std::sqrt( ( differenceRightAscensionEuropaGanymede * std::cos( currentDeclinationEuropa ) )
                                                                  * ( differenceRightAscensionEuropaGanymede * std::cos( currentDeclinationEuropa ) )
                                                     + differenceDeclinationEuropaGanymede * differenceDeclinationEuropaGanymede );

        Eigen::VectorXd currentApparentDistancesVector = Eigen::VectorXd::Zero( 3 );
        currentApparentDistancesVector[ 0 ] = currentApparentDistanceIoEuropa * 180.0 / tudat::mathematical_constants::PI * 3600.0;
        currentApparentDistancesVector[ 1 ] = currentApparentDistanceIoGanymede * 180.0 / tudat::mathematical_constants::PI * 3600.0;
        currentApparentDistancesVector[ 2 ] = currentApparentDistanceEuropaGanymede * 180.0 / tudat::mathematical_constants::PI * 3600.0;

        apparentDistancesHistory[ itr->first ] = currentApparentDistancesVector;

        Eigen::VectorXd currentAngularPositionsVector = Eigen::VectorXd::Zero( 6 );
        currentAngularPositionsVector = ( Eigen::Vector6d( ) << currentRightAscensionIo, currentDeclinationIo,
                                          currentRightAscensionEuropa, currentDeclinationEuropa,
                                          currentRightAscensionGanymede, currentDeclinationGanymede ).finished( );
        angularPositionsHistory[ itr->first ] = currentAngularPositionsVector;

    }


    input_output::writeDataMapToTextFile( apparentDistancesHistory,
                                          "apparentDistancesHistory_SecondTestCase.dat",
                                          "C:/Users/chamb/Documents/PhD/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    Eigen::Vector6d initialStateJupiter = newSystemInitialState.segment( 6, 6 );
    Eigen::Vector6d initialStateIo = newSystemInitialState.segment( 12, 6 );
    Eigen::Vector6d initialStateEuropa = newSystemInitialState.segment( 18, 6 );
    Eigen::Vector6d initialStateGanymede = newSystemInitialState.segment( 24, 6 );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Earth", initialStateJupiter, "SSB" ) );
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Jupiter", initialStateJupiter, "SSB" ) );
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Io", initialStateIo, "Jupiter" ) );
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Europa", initialStateEuropa, "Jupiter" ) );
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Ganymede", initialStateGanymede, "Jupiter" ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );


    double estimatedCentralInstant = 1.341664005849638e6; // from Matlab


    // Create light-time correction settings
    std::vector< std::string > perturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                perturbingBodies ) );

    // Create observation settings
    std::shared_ptr< observation_models::ObservationSettings > observableSettings = std::make_shared< observation_models::ObservationSettings >
            ( observation_models::apparent_distance, lightTimeCorrectionSettings, std::make_shared< observation_models::ConstantObservationBiasSettings >(
                  ( Eigen::Vector1d( ) << 0.0 ).finished( ), true ) );

    observation_models::LinkEnds linkEnds;
    linkEnds[ observation_models::receiver ] = std::make_pair( "Earth" , ""  );
    linkEnds[ observation_models::transmitter ] = std::make_pair( "Io" , ""  );
    linkEnds[ observation_models::transmitter2 ] = std::make_pair( "Europa" , ""  );

    // Create apparent distance observation model.
    std::shared_ptr< observation_models::ObservationModel< 1, double, double > > apparentDistanceModel =
           observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                linkEnds, observableSettings, bodyMap );

    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector1d observationFromReceptionTime = 3600.0 * 180.0 / mathematical_constants::PI * apparentDistanceModel->computeObservationsWithLinkEndData( estimatedCentralInstant, observation_models::receiver, linkEndTimes, linkEndStates );
    std::cout << "apparent distance at central instant from observation model: " << observationFromReceptionTime << "\n\n";
    std::cout << "link end times: " << linkEndTimes[ 0 ] << " & " << linkEndTimes[ 1 ] << " & " << linkEndTimes[ 2 ] << "\n\n";

    std::map< double, Eigen::VectorXd > lightTimeCorrectedApparentDistances;
    for ( double i = estimatedCentralInstant - 20.0 * 60.0 ; i <= estimatedCentralInstant + 30.0 * 60.0 ; i += 30.0 )
    {
        lightTimeCorrectedApparentDistances[ i ] = 3600.0 * 180.0 / mathematical_constants::PI *
                apparentDistanceModel->computeObservations( i, observation_models::receiver );
    }


    input_output::writeDataMapToTextFile( lightTimeCorrectedApparentDistances,
                                          "lightTimeCorrectedApparentDistancesHistory.dat",
                                          "C:/Users/chamb/Documents/PhD/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    // Define more complete list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > updatedDependentVariablesList;
    updatedDependentVariablesList = dependentVariablesList;

    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Io", "Sun" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Io", "Jupiter" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Io", "Europa" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Io", "Ganymede" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Europa", "Sun" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Europa", "Jupiter" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Europa", "Ganymede" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_velocity_dependent_variable, "Io", "Earth" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_velocity_dependent_variable, "Europa", "Earth" ) );
    updatedDependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                                 relative_position_dependent_variable, "Earth", "Sun" ) );

    updatedDependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Io" ) );
    updatedDependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Europa" ) );
    updatedDependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, "Earth" ) );

    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Io", "Jupiter" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Earth", "Jupiter" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Io", "Europa", "Jupiter" ) );

    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Europa", "Jupiter" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Earth", "Jupiter" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Europa", "Io", "Jupiter" ) );

    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Earth", "Sun" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Io", "Sun" ) );
    updatedDependentVariablesList.push_back(
                std::make_shared< propagators::TotalAccelerationPartialWrtStateSaveSettings >(
                    "Earth", "Europa", "Sun" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > updatedDependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( updatedDependentVariablesList );



    // Redefine propagator settings.
    newPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, newSystemInitialState,
                std::make_shared< propagators::PropagationTimeTerminationSettings > ( estimatedCentralInstant, true ),
                cowell, updatedDependentVariablesToSave );

    // Create simulation object and propagate dynamics.
    SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                bodyMap, newIntegratorSettings, newPropagatorSettings, parametersToEstimate, true,
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true );


    // Retrieve dependent variables history.
    std::map< double, Eigen::VectorXd > dependentVariablesHistoryForInterface = variationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );
    // Create dependent variables interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistoryForInterface ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistoryForInterface ), 4 );

    std::shared_ptr< propagators::SingleArcDependentVariablesInterface > dependentVariablesInterface = std::make_shared< propagators::SingleArcDependentVariablesInterface >(
                dependentVariablesInterpolator, updatedDependentVariablesToSave );

    std::shared_ptr< observation_partials::MutualApproximationScaling > mutualApproximationScaling =
            std::make_shared< observation_partials::MutualApproximationScaling >( dependentVariablesInterface );




    Eigen::VectorXd centralInstantRelativePositionIoEarth = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getDependentVariableHistory( ).rbegin( )->second.segment( 0, 3 );
    Eigen::VectorXd centralInstantRelativePositionEuropaEarth = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getDependentVariableHistory( ).rbegin( )->second.segment( 3, 3 );

    Eigen::VectorXd centralInstantDependentVariables = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getDependentVariableHistory( ).rbegin( )->second;

    Eigen::Matrix< double, 3, 1 > centralInstantSphericalStateIo =
            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                centralInstantRelativePositionIoEarth.segment( 0, 3 ) ).template cast< double >( );
    double centralInstantRightAscensionIo = centralInstantSphericalStateIo.z( );
    double centralInstantDeclinationIo = mathematical_constants::PI / 2.0 - centralInstantSphericalStateIo.y( );

    Eigen::Matrix< double, 3, 1 > centralInstantSphericalStateEuropa =
            tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                centralInstantRelativePositionEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
    double centralInstantRightAscensionEuropa = centralInstantSphericalStateEuropa.z( );
    double centralInstantDeclinationEuropa = mathematical_constants::PI / 2.0 - centralInstantSphericalStateEuropa.y( );

    double centralInstantDifferenceRightAscensionIoEuropa = centralInstantRightAscensionEuropa - centralInstantRightAscensionIo;
    double centralInstantDifferenceDeclinationIoEuropa = centralInstantDeclinationEuropa - centralInstantDeclinationIo;
    double centralInstantApparentDistanceIoEuropa =
            std::sqrt( ( centralInstantDifferenceRightAscensionIoEuropa * std::cos( centralInstantDeclinationIo ) )
                       * ( centralInstantDifferenceRightAscensionIoEuropa * std::cos( centralInstantDeclinationIo ) )
                       + centralInstantDifferenceDeclinationIoEuropa * centralInstantDifferenceDeclinationIoEuropa );

    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
    std::cout << "apparent distance between Io and Europa at estimated central instant: " <<
                 centralInstantApparentDistanceIoEuropa * 180.0 / mathematical_constants::PI * 3600.0 << "\n\n";

    std::map< double, Eigen::MatrixXd > stateTransitionResult =
            variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 0 );

    Eigen::VectorXd centralInstantState = variationalEquationsSimulator.getDynamicsSimulator( )
            ->getEquationsOfMotionNumericalSolution( ).rbegin( )->second;

    Eigen::MatrixXd centralInstantStateTransitionMatrix = stateTransitionResult.rbegin( )->second;
    Eigen::Vector3d centralInstantVectorFromIoToEuropa = centralInstantState.segment( 18, 3 ) - centralInstantState.segment( 12, 3 );
    Eigen::Vector3d centralInstantUnitVectorFromIoToEuropa = centralInstantVectorFromIoToEuropa / centralInstantVectorFromIoToEuropa.norm( );
    double centralInstantDistanceIoEuropa = centralInstantVectorFromIoToEuropa.norm( );

    Eigen::Vector3d variationPositionCentralInstant = centralInstantDistanceIoEuropa * 0.005
            * centralInstantUnitVectorFromIoToEuropa;


    for ( unsigned int testCase = 0 ; testCase < 9 ; testCase++ )
    {

        Eigen::VectorXd totalVariationVectorCentralInstant = Eigen::VectorXd::Zero( 6 * 5 );
        std::string nameOutputFile;

        if ( testCase == 0 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << centralInstantState[ 12 ] * 0.0001, 0.0 , 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 12, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationXio.dat";
        }
        else if ( testCase == 1 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, centralInstantState[ 13 ] * 0.0001, 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 12, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationYio.dat";
        }
        else if ( testCase == 2 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, 0.0, centralInstantState[ 14 ] * 0.0001 ).finished( );
            totalVariationVectorCentralInstant.segment( 12, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationZio.dat";
        }
        else if ( testCase == 3 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << centralInstantState[ 18 ] * 0.0001, 0.0, 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 18, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationXeuropa.dat";
        }
        else if ( testCase == 4 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, centralInstantState[ 19 ] * 0.0001, 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 18, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationYeuropa.dat";
        }
        else if ( testCase == 5 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, 0.0, centralInstantState[ 20 ] * 0.0001 ).finished( );
            totalVariationVectorCentralInstant.segment( 18, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationZeuropa.dat";
        }
        else if ( testCase == 6 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << centralInstantState[ 0 ] * 0.0001, 0.0, 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 0, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationXearth.dat";
        }
        else if ( testCase == 7 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, centralInstantState[ 1 ] * 0.0001, 0.0 ).finished( );
            totalVariationVectorCentralInstant.segment( 0, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationYearth.dat";
        }
        else if ( testCase == 8 )
        {
            variationPositionCentralInstant = ( Eigen::Vector3d( ) << 0.0, 0.0, centralInstantState[ 2 ] * 0.0001 ).finished( );
            totalVariationVectorCentralInstant.segment( 0, 3 ) = variationPositionCentralInstant;
            nameOutputFile = "apparentDistancesHistory_variationZearth.dat";
        }

        Eigen::VectorXd variationVectorInitialEpoch = centralInstantStateTransitionMatrix.inverse( ) * totalVariationVectorCentralInstant;

        std::cout << "variation vector initial epoch: " << variationVectorInitialEpoch.transpose( ) << "\n\n";


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE PERTURBED INITIAL STATE                 /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        Eigen::VectorXd updatedSystemInitialState = newSystemInitialState + variationVectorInitialEpoch;

        // Define propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > updatedPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > > (
                    centralBodies, accelerationModelMap, bodiesToPropagate, updatedSystemInitialState,
                    newCustomTerminationSettings, cowell, dependentVariablesToSave );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > updatedDynamicsSimulator( bodyMap, newIntegratorSettings, updatedPropagatorSettings, true, false, true );
        std::map< double, Eigen::VectorXd > updatedDependentVariablesHistory = updatedDynamicsSimulator.getDependentVariableHistory( );


        // Create map with apparent relative distances.
        std::map< double, Eigen::VectorXd > updatedApparentDistancesHistory;

        for ( std::map< double, Eigen::VectorXd >::iterator itr = updatedDependentVariablesHistory.begin( ) ;
              itr != updatedDependentVariablesHistory.end( ) ; itr++ )
        {
            Eigen::VectorXd relativePositionIoEarth = itr->second.segment( 0, 3 );
            Eigen::VectorXd relativePositionEuropaEarth = itr->second.segment( 3, 3 );
            Eigen::VectorXd relativePositionGanymedeEarth = itr->second.segment( 6, 3 );

            Eigen::Matrix< double, 3, 1 > sphericalStateIo = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                        relativePositionIoEarth.segment( 0, 3 ) ).template cast< double >( );
            double currentRightAscensionIo = sphericalStateIo.z( );
            double currentDeclinationIo = mathematical_constants::PI / 2.0 - sphericalStateIo.y( );

            Eigen::Matrix< double, 3, 1 > sphericalStateEuropa = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                        relativePositionEuropaEarth.segment( 0, 3 ) ).template cast< double >( );
            double currentRightAscensionEuropa = sphericalStateEuropa.z( );
            double currentDeclinationEuropa = mathematical_constants::PI / 2.0 - sphericalStateEuropa.y( );

            Eigen::Matrix< double, 3, 1 > sphericalStateGanymede = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                        relativePositionGanymedeEarth.segment( 0, 3 ) ).template cast< double >( );
            double currentRightAscensionGanymede = sphericalStateGanymede.z( );
            double currentDeclinationGanymede = mathematical_constants::PI / 2.0 - sphericalStateGanymede.y( );

            double differenceRightAscensionIoEuropa = currentRightAscensionEuropa - currentRightAscensionIo;
            double differenceDeclinationIoEuropa = currentDeclinationEuropa - currentDeclinationIo;
            double currentApparentDistanceIoEuropa = std::sqrt( ( differenceRightAscensionIoEuropa * std::cos( currentDeclinationIo ) )
                                                                * ( differenceRightAscensionIoEuropa * std::cos( currentDeclinationIo ) )
                                                         + differenceDeclinationIoEuropa * differenceDeclinationIoEuropa );

            double differenceRightAscensionIoGanymede = currentRightAscensionGanymede - currentRightAscensionIo;
            double differenceDeclinationIoGanymede = currentDeclinationGanymede - currentDeclinationIo;
            double currentApparentDistanceIoGanymede = std::sqrt( ( differenceRightAscensionIoGanymede * std::cos( currentDeclinationIo ) )
                                                                  * ( differenceRightAscensionIoGanymede * std::cos( currentDeclinationIo ) )
                                                         + differenceDeclinationIoGanymede * differenceDeclinationIoGanymede );

            double differenceRightAscensionEuropaGanymede = currentRightAscensionGanymede - currentRightAscensionEuropa;
            double differenceDeclinationEuropaGanymede = currentDeclinationGanymede - currentDeclinationEuropa;
            double currentApparentDistanceEuropaGanymede = std::sqrt( ( differenceRightAscensionEuropaGanymede * std::cos( currentDeclinationEuropa ) )
                                                                      * ( differenceRightAscensionEuropaGanymede * std::cos( currentDeclinationEuropa ) )
                                                         + differenceDeclinationEuropaGanymede * differenceDeclinationEuropaGanymede );

            Eigen::VectorXd currentApparentDistancesVector = Eigen::VectorXd::Zero( 3 );
            currentApparentDistancesVector[ 0 ] = currentApparentDistanceIoEuropa * 180.0 / tudat::mathematical_constants::PI * 3600.0;
            currentApparentDistancesVector[ 1 ] = currentApparentDistanceIoGanymede * 180.0 / tudat::mathematical_constants::PI * 3600.0;
            currentApparentDistancesVector[ 2 ] = currentApparentDistanceEuropaGanymede * 180.0 / tudat::mathematical_constants::PI * 3600.0;

            updatedApparentDistancesHistory[ itr->first ] = currentApparentDistancesVector;

        }


        input_output::writeDataMapToTextFile( updatedApparentDistancesHistory,
                                              nameOutputFile,
                                              "C:/Users/chamb/Documents/PhD/",
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

    }


    Eigen::Vector3d centralInstantEarthToIo = centralInstantDependentVariables.segment( 0, 3 );
    Eigen::Vector3d centralInstantEarthToEuropa = centralInstantDependentVariables.segment( 3, 3 );

    Eigen::Vector3d centralInstantRelativeVelocityIoWrtEarth = centralInstantDependentVariables.segment( 33, 3 );
    Eigen::Vector3d centralInstantRelativeVelocityEuropaWrtEarth = centralInstantDependentVariables.segment( 36, 3 );


    double numericalEstimationPartialT0wrtXio = -0.865298238582909 / ( centralInstantState[ 12 ] * 0.0001 );
    double numericalEstimationPartialT0wrtYio = -2.822315546218306 / ( centralInstantState[ 13 ] * 0.0001 );
    double numericalEstimationPartialT0wrtZio = -0.540046634851024 / ( centralInstantState[ 14 ] * 0.0001 );

    double numericalEstimationPartialT0wrtXeuropa = 2.921673136996105 / ( centralInstantState[ 18 ] * 0.0001 );
    double numericalEstimationPartialT0wrtYeuropa = 1.120536909671500 / ( centralInstantState[ 19 ] * 0.0001 );
    double numericalEstimationPartialT0wrtZeuropa = 0.178772800834849 / ( centralInstantState[ 20 ] * 0.0001 );

    double numericalEstimationPartialT0wrtXearth = -0.207597788888961 / ( centralInstantState[ 0 ] * 0.0001 );
    double numericalEstimationPartialT0wrtYearth = -0.728156324941665 / ( centralInstantState[ 1 ] * 0.0001 );
    double numericalEstimationPartialT0wrtZearth = -0.111029546242207 / ( centralInstantState[ 2 ] * 0.0001 );


    // Update mutual approximation scaling.
    linkEndStates.clear( );
    Eigen::Vector6d currentStateForLinkEnd;
    currentStateForLinkEnd.segment( 0, 3 ) = centralInstantEarthToIo;
    currentStateForLinkEnd.segment( 3, 3 ) = centralInstantRelativeVelocityIoWrtEarth;
    linkEndStates.push_back( currentStateForLinkEnd );
    currentStateForLinkEnd.segment( 0, 3 ) = centralInstantEarthToEuropa;
    currentStateForLinkEnd.segment( 3, 3 ) = centralInstantRelativeVelocityEuropaWrtEarth;
    linkEndStates.push_back( currentStateForLinkEnd );
    linkEndStates.push_back( Eigen::Vector6d::Zero( ) );

    std::vector< double > times;
    times.push_back( estimatedCentralInstant );
    times.push_back( estimatedCentralInstant );
    times.push_back( estimatedCentralInstant );

    Eigen::VectorXd currentObservation;
    mutualApproximationScaling->update( linkEndStates, times, observation_models::receiver, linkEnds, currentObservation );
    Eigen::Matrix< double, 1, 3 > partialCentralInstantWrtFirstTransmitterPosition = mutualApproximationScaling->getScalingFactor( observation_models::transmitter );
    Eigen::Matrix< double, 1, 3 > partialCentralInstantWrtSecondTransmitterPosition = mutualApproximationScaling->getScalingFactor( observation_models::transmitter2 );
    Eigen::Matrix< double, 1, 3 > partialCentralInstantWrtReceiverPosition = mutualApproximationScaling->getScalingFactor( observation_models::receiver );

    std::cout << "partials central instant w.r.t. Io position from getScalingFactor function: " << partialCentralInstantWrtFirstTransmitterPosition << "\n\n";
    std::cout << "numerical partials central instant w.r.t. Io position: " << numericalEstimationPartialT0wrtXio << "  " <<
              numericalEstimationPartialT0wrtYio << "  " << numericalEstimationPartialT0wrtZio << "\n\n";

    std::cout << "partials central instant w.r.t. Europa position from getScalingFactor function: " << partialCentralInstantWrtSecondTransmitterPosition << "\n\n";
    std::cout << "numerical partials central instant w.r.t. Europa position: " << numericalEstimationPartialT0wrtXeuropa << "  " <<
              numericalEstimationPartialT0wrtYeuropa << "  " << numericalEstimationPartialT0wrtZeuropa << "\n\n";

    std::cout << "partials central instant w.r.t. Earth position from getScalingFactor function: " << partialCentralInstantWrtReceiverPosition << "\n\n";
    std::cout << "numerical partials central instant w.r.t. Earth position: " << numericalEstimationPartialT0wrtXearth << "  " <<
              numericalEstimationPartialT0wrtYearth << "  " << numericalEstimationPartialT0wrtZearth << "\n\n";

    std::cout << "relative difference partial t0 wrt x Io: " <<
                 ( partialCentralInstantWrtFirstTransmitterPosition[ 0 ] - numericalEstimationPartialT0wrtXio ) / numericalEstimationPartialT0wrtXio << "\n\n";
    std::cout << "relative difference partial t0 wrt y Io: " <<
                 ( partialCentralInstantWrtFirstTransmitterPosition[ 1 ] - numericalEstimationPartialT0wrtYio ) / numericalEstimationPartialT0wrtYio << "\n\n";
    std::cout << "relative difference partial t0 wrt z Io: " <<
                 ( partialCentralInstantWrtFirstTransmitterPosition[ 2 ] - numericalEstimationPartialT0wrtZio ) / numericalEstimationPartialT0wrtZio << "\n\n";


    std::cout << "relative difference partial t0 wrt x Europa: " <<
                 ( partialCentralInstantWrtSecondTransmitterPosition[ 0 ] - numericalEstimationPartialT0wrtXeuropa ) / numericalEstimationPartialT0wrtXeuropa << "\n\n";
    std::cout << "relative difference partial t0 wrt y Europa: " <<
                 ( partialCentralInstantWrtSecondTransmitterPosition[ 1 ] - numericalEstimationPartialT0wrtYeuropa ) / numericalEstimationPartialT0wrtYeuropa << "\n\n";
    std::cout << "relative difference partial t0 wrt z Europa: " <<
                 ( partialCentralInstantWrtSecondTransmitterPosition[ 2 ] - numericalEstimationPartialT0wrtZeuropa ) / numericalEstimationPartialT0wrtZeuropa << "\n\n";


    std::cout << "relative difference partial t0 wrt x Earth: " <<
                 ( partialCentralInstantWrtReceiverPosition[ 0 ] - numericalEstimationPartialT0wrtXearth ) / numericalEstimationPartialT0wrtXearth << "\n\n";
    std::cout << "relative difference partial t0 wrt y Earth: " <<
                 ( partialCentralInstantWrtReceiverPosition[ 1 ] - numericalEstimationPartialT0wrtYearth ) / numericalEstimationPartialT0wrtYearth << "\n\n";
    std::cout << "relative difference partial t0 wrt z Earth: " <<
                 ( partialCentralInstantWrtReceiverPosition[ 2 ] - numericalEstimationPartialT0wrtZearth ) / numericalEstimationPartialT0wrtZearth << "\n\n";


}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





