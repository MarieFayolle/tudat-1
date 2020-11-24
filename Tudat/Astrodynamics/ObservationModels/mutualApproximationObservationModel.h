/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALAPPROXIMATIONOBSERVATIONMODEL_H
#define TUDAT_MUTUALAPPROXIMATIONOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/apparentDistanceObservationModel.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"
#include "Tudat/Mathematics/RootFinders/createRootFinder.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/mutualApproximationPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/observationViabilityCalculator.h"

namespace tudat
{

namespace observation_models
{


//! Polynomial function for the root-finder.
struct PolynomialFunction : public basic_mathematics::BasicFunction< double, double >
{

    //! Create a function that returns the time of flight associated with the current trajectory.
    PolynomialFunction( const Eigen::VectorXd polynomialCoefficients,
                        const unsigned int orderPolynomialFunction ):
        polynomialCoefficients_( polynomialCoefficients ),
        orderPolynomialFunction_( orderPolynomialFunction )
    {
        if ( orderPolynomialFunction_ != ( polynomialCoefficients_.size( ) - 1 ) )
        {
            throw std::runtime_error( "Error when computing the roots of a polynomial function, number of polynomial coefficients unconsistent with the order"
                                      "for the polynomial function." );
        }
    }

    //! Evaluate the difference between the current and required time of flight values, from the current value of the free coefficient.
    double evaluate( const double inputValue )
    {
        double outputValue = 0.0;
        for ( unsigned int order = 0 ; order <= orderPolynomialFunction_ ; order++ )
        {
            outputValue += polynomialCoefficients_[ order ] * std::pow( inputValue, order );
        }
        return outputValue;
    }

    //! Derivative function (not provided).
    double computeDerivative( const unsigned int order, const double inputValue )
    {
        // Return zero if the order of the derivative is larger than the order of the polynomial function.
        if ( order > orderPolynomialFunction_ )
        {
            return 0.0;
        }
        else
        {
            double outputValue = 0.0;
            for ( unsigned int i = orderPolynomialFunction_ ; i >= order ; i-- )
            {
                double derivativeCoefficient = polynomialCoefficients_[ i ];
                for ( unsigned int j = 0 ; j < order ; j++ )
                {
                    derivativeCoefficient *= ( i - j );
                }
                outputValue += derivativeCoefficient * std::pow( inputValue, i - order );
            }
            return outputValue;
        }
    }

    //! Integral function (not provided).
    double computeDefiniteIntegral( unsigned int order, double lowerBound, double upperbound )
    {
        throw std::runtime_error( "The rootfinder for polynomial function should not evaluate integrals!" );
    }

    //! Get the expected true location of the root (not implemented).
    double getTrueRootLocation( ) { return TUDAT_NAN; }

    //! Get the accuracy of the true location of the root (not implemented).
    double getTrueRootAccuracy( ) { return TUDAT_NAN; }

    //! Get a reasonable initial guess of the root location (not implemented).
    double getInitialGuess( ) { return TUDAT_NAN; }

    //! Get a reasonable lower boundary for the root location (not implemented).
    double getLowerBound( ) { return TUDAT_NAN; }

    //! Get a reasonable upper boundary for the root location (not implemented).
    double getUpperBound( ) { return TUDAT_NAN; }

protected:

private:

    Eigen::VectorXd polynomialCoefficients_;
    unsigned int orderPolynomialFunction_;
};

//! Class for simulating mutual approximation (central instant) observables.
/*!
 *  Class for simulating mutual approximation (central instant), using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class MutualApproximationObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     * Constructor
     * \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the first transmitter and the receiver.
     * \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the second transmitter and the receiver.
     * \param frequencyApparentDistanceObservations Time interval between two successive apparent distance measurements.
     * \param toleranceWrtExpectedCentralInstant Tolerance w.r.t. to the estimated central instant (time interval over which the central instant is looked for).
     * \param upperLimitImpactParameter Upper limit value for the impact parameter, beyond which there is no close encounter.
     * \param orderPolynomialFitting Order of the fitting polynomial to derive the central instant.
     * \param checkExistenceMutualApproximation Boolean denoting whether the existence of a mutual approximation at the given time must be verified.
     * \param checkObservableDuplicates Boolean denoting whether the duplicated central instant observables should be removed.
     * \param rootFinderSettings Root finder settings used to derive the central instant from the apparent distance observations.
     * \param mutualApproximationBiasCalculator Object for calculating system-dependent errors in the central instant
     *  observable, i.e. deviations from the physically ideal central instant observable between reference points (default none).
     * \param apparentDistancesBiasCalculator Object for calculating system-dependent errors in the apparent distance
     *  observables, i.e. deviations from the physically ideal apparent distance observables between reference points (default none).
     */
    MutualApproximationObservationModel(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter,
            const double frequencyApparentDistanceObservations,
            const double toleranceWrtExpectedCentralInstant,
            const double upperLimitImpactParameter,
            const int orderPolynomialFitting = 4,
            const bool checkExistenceMutualApproximation = true,
            const bool checkObservableDuplicates = true,
            const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
            const std::shared_ptr< ObservationBias< 1 > > mutualApproximationBiasCalculator = nullptr,
            const std::shared_ptr< ObservationBias< 1 > > apparentDistancesBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( mutual_approximation, mutualApproximationBiasCalculator ),
        lightTimeCalculatorFirstTransmitter_( lightTimeCalculatorFirstTransmitter ),
        lightTimeCalculatorSecondTransmitter_( lightTimeCalculatorSecondTransmitter ),
        frequencyApparentDistanceObservations_( frequencyApparentDistanceObservations ),
        toleranceWrtExpectedCentralInstant_( toleranceWrtExpectedCentralInstant ),
        upperLimitImpactParameter_( upperLimitImpactParameter ),
        orderPolynomialFitting_( orderPolynomialFitting ),
        checkExistenceMutualApproximation_( checkExistenceMutualApproximation ),
        checkObservableDuplicates_( checkObservableDuplicates ),
        rootFinderSettings_( rootFinderSettings )/*,
        viabilityCalculatorsApparentDistances_( viabilityCalculatorsApparentDistances )*/
    {
        // Create required observation model to compute the intermediate apparent distance measurements.
        apparentDistanceObservationModel_ =
                std::make_shared< ApparentDistanceObservationModel< ObservationScalarType, TimeType > >
                ( lightTimeCalculatorFirstTransmitter_, lightTimeCalculatorSecondTransmitter_, apparentDistancesBiasCalculator );

        // Initialise empty vector of viability calculator for intermediate apparent distance observations.
        viabilityCalculatorsApparentDistances_ = std::vector< std::shared_ptr< ObservationViabilityCalculator > >( );

    }

    //! Destructor
    ~MutualApproximationObservationModel( ){ }

    //! Function to compute ideal central instants for mutual approximations.
    /*!
     *  This function compute ideal central instants for mutual approximations.
     *  The time argument should always be the reception time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which the mutual approximation is expected to occur.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value). Should always be receiver.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated central instant observable values.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        std::cout.precision( 20 );

        double currentTime = time;

        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating mutual approximation observation, link end is not receiver." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState;

        ObservationScalarType estimatedCentralInstant = TUDAT_NAN;
        bool iterativeProcess = true;
        unsigned int iteration = 1;
        while ( iterativeProcess && ( iterativeProcess <= 2 ) ){

            if ( iteration == 2 )
            {
                iterativeProcess = false;
            }

        // Compute intermediate apparent distance observables.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations;

        double observationStartingTime = currentTime - toleranceWrtExpectedCentralInstant_;
        double observationEndingTime = currentTime + toleranceWrtExpectedCentralInstant_;
//        double currentObservationTime = observationStartingTime;

        for ( TimeType currentObservationTime = observationStartingTime ;
              currentObservationTime <= observationEndingTime ;
              currentObservationTime += frequencyApparentDistanceObservations_ )
        {
//            Eigen::Vector1d apparentDistanceObservation = apparentDistanceObservationModel_->computeObservations( currentObservationTime, receiver );
//            std::cout << "current observation time and apparent distance: " << currentObservationTime
//                      << " & " << apparentDistanceObservation * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceObservation =
                    apparentDistanceObservationModel_->computeObservationsWithLinkEndData( currentObservationTime, receiver, vectorOfTimes, vectorOfStates );

            bool observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, viabilityCalculatorsApparentDistances_ );
            if ( observationFeasible )
            {
                apparentDistanceObservations[ currentObservationTime ] = apparentDistanceObservation;
            }
        }


//        /// TO BE REMOVED EVENTUALLY
//        bool viableCentralInstant = isCentralInstantObservableViable( time, upperLimitImpactParameter_, apparentDistanceObservations );
//        std::cout << "is mutual approximation observation viable " << viableCentralInstant << "\n\n";

        /*ObservationScalarType*/ estimatedCentralInstant = TUDAT_NAN;
        if ( ( !checkExistenceMutualApproximation_ ) || ( isCentralInstantObservableViable( currentTime, upperLimitImpactParameter_, apparentDistanceObservations ) ) )
        {

    //        std::cout << "size apparent distance observations map: " << apparentDistanceObservations.size( ) << "\n\n";


            // Polynomial fitting.
            std::vector< double > polynomialPowers; // = { 0.0, 1.0, 2.0, 3.0, 4.0 };
            for ( unsigned int i = 0 ; i <= orderPolynomialFitting_ ; i++ )
            {
                polynomialPowers.push_back( i );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistances;
            apparentDistances.resize( apparentDistanceObservations.size( ) );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistanceTimes;
            apparentDistanceTimes.resize( apparentDistanceObservations.size( ) );

            double initialisationCounter = - ( int ) ( apparentDistanceObservations.size( ) / 2.0 );

    //        std::cout << "TEST: " << - ( int ) ( apparentDistanceObservations.size( ) / 2.0 ) << "\n\n";

            typename std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > >::iterator itr = apparentDistanceObservations.begin( );
            for ( unsigned int i = 0 ; i < apparentDistanceObservations.size( ) ; i++ )
            {
                apparentDistanceTimes( i, 0 ) = i + initialisationCounter; //itr->first / 3600.0;
                apparentDistances( i, 0 ) =  itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI;
                itr++;
            }

    //        std::cout << "apparent distance times: " << apparentDistanceTimes.transpose( ) << "\n\n";
    //        std::cout << "apparent distances: " << apparentDistances.transpose( ) << "\n\n";


            Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceTimes, apparentDistances, polynomialPowers );

            Eigen::VectorXd derivativePolynomialCoefficients; derivativePolynomialCoefficients.resize( polynomialCoefficients.size( ) - 1 );
            for ( unsigned int i = 1 ; i <= orderPolynomialFitting_ ; i++ )
            {
                derivativePolynomialCoefficients[ i - 1 ] = i * polynomialCoefficients[ i ];
            }

//            std::cout << "polynomial coefficients: " << polynomialCoefficients.transpose( ) << "\n\n";

    //        std::cout << "derivative polynomial coefficients: " << derivativePolynomialCoefficients.transpose( ) << "\n\n";

    //        double counter = initialisationCounter;
    //        for ( typename std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > >::iterator itr = apparentDistanceObservations.begin( ) ;
    //              itr != apparentDistanceObservations.end( ) ; itr++ )
    //        {
    //            std::cout << "true apparent distance: " << itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
    //            std::cout << "polynomial approximation: " <<
    //                         polynomialCoefficients[ 0 ] + polynomialCoefficients[ 1 ] * counter
    //                    + polynomialCoefficients[ 2 ] * ( counter * counter )
    //                    + polynomialCoefficients[ 3 ] * ( counter * counter * counter )
    //                    + polynomialCoefficients[ 4 ] * ( counter * counter * counter * counter ) << "\n\n";
    //            std::cout << "residuals: " << itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI
    //                         - ( polynomialCoefficients[ 0 ] + polynomialCoefficients[ 1 ] * counter
    //                         + polynomialCoefficients[ 2 ] * ( counter * counter )
    //                         + polynomialCoefficients[ 3 ] * ( counter * counter * counter )
    //                         + polynomialCoefficients[ 4 ] * ( counter * counter * counter * counter ) ) << "\n\n";
    //            counter += 1;
    //        }

            double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
            double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
            double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];

            double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
            double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

            double beta = Q * Q * Q + R * R;



    //        ObservationScalarType estimatedCentralInstant;
            if ( beta < 0 )
            {
                double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
                double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
                double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
                double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "t2: " << t2 << "\n\n";
//                std::cout << "t3: " << t3 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t2 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t3 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
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
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
            }

//            std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";

            // Use of root finder given as input.
            std::shared_ptr< basic_mathematics::Function< double, double > > derivativePolynomialFittingFunction =
                    std::make_shared< PolynomialFunction >( derivativePolynomialCoefficients, orderPolynomialFitting_ - 1 );

            // Create root finder from root finder settings.
            double estimatedCentralInstant2 = TUDAT_NAN;
            try
            {
                std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder
                        = root_finders::createRootFinder( rootFinderSettings_, initialisationCounter, - initialisationCounter, 0.0 );
//                estimatedCentralInstant2 = observationStartingTime
//                    + ( rootFinder->execute( derivativePolynomialFittingFunction, 0.0 ) - initialisationCounter ) * frequencyApparentDistanceObservations_;

//                std::cout << "estimated central instant root finder: " << estimatedCentralInstant2 << "\n\n";
//                std::cout << "difference in central instant between analytical and root finder solutions: "
//                          << estimatedCentralInstant - estimatedCentralInstant2 << "\n\n";

        //        estimatedCentralInstant2 = rootFinder->execute( derivativePolynomialFittingFunction, 0.0 );
        //        std::cout << "TEST: " << derivativePolynomialCoefficients[ 0 ]
        //                     + derivativePolynomialCoefficients[ 1 ] * estimatedCentralInstant2
        //                + derivativePolynomialCoefficients[ 2 ] * ( estimatedCentralInstant2 * estimatedCentralInstant2 )
        //                + derivativePolynomialCoefficients[ 3 ] * ( estimatedCentralInstant2 * estimatedCentralInstant2 * estimatedCentralInstant2 ) << "\n\n";
        //        double estimatedCentralInstant3 = ( estimatedCentralInstant - observationStartingTime ) / frequencyApparentDistanceObservations_
        //                + initialisationCounter;
        //        std::cout << "TEST2: " << derivativePolynomialCoefficients[ 0 ]
        //                     + derivativePolynomialCoefficients[ 1 ] * estimatedCentralInstant3
        //                + derivativePolynomialCoefficients[ 2 ] * ( estimatedCentralInstant3 * estimatedCentralInstant3 )
        //                + derivativePolynomialCoefficients[ 3 ] * ( estimatedCentralInstant3 * estimatedCentralInstant3 * estimatedCentralInstant3 ) << "\n\n";

                // Compute light-time and receiver/transmitter states.
                ObservationScalarType lightTimeFirstTransmitter = lightTimeCalculatorFirstTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, firstTransmitterState, estimatedCentralInstant, true );

                ObservationScalarType lightTimeSecondTransmitter = lightTimeCalculatorSecondTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, secondTransmitterState, estimatedCentralInstant, true );


                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndStates.clear( );

                linkEndStates.push_back( firstTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( secondTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( receiverState.template cast< double >( ) );

                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeFirstTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeSecondTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant ) );

                if ( checkObservableDuplicates_ )
                {
                    if ( !isMutualApproximationAlreadyDetected( estimatedCentralInstant ) )
                    {
                        if ( !iterativeProcess )
                        {
                        centralInstantsOfDetectedMutualApproximations_.push_back( estimatedCentralInstant );
                        }
                    }
                    else
                    {
                        estimatedCentralInstant = TUDAT_NAN;

                        // Set link end times and states.
                        linkEndTimes.clear( );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );

                        iterativeProcess = false;
                    }
                }
            }
            catch ( std::runtime_error& error )
            {
                std::cerr << "Warning, no mutual approximation is found around the requested observation time t = " << std::to_string( currentTime ) <<
                             ". The root finder has failed to find a minimum in the apparent distances history. Returns NAN as observable value." << std::endl;
                estimatedCentralInstant = TUDAT_NAN;

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );

                linkEndStates.clear( );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

                iterativeProcess = false;
            }

        }
        else
        {
            // Set link end times and states.
            linkEndTimes.clear( );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );

            linkEndStates.clear( );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

            iterativeProcess = false;
        }


        iteration++;
        currentTime = estimatedCentralInstant;

    }



        /// Compute slightly modified observable for different formulation of the variational equations.

        // Compute spherical relative position for first transmitter wrt receiver.
        Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinatesFirstTransmitter =
                coordinate_conversions::convertCartesianToSpherical< ObservationScalarType >(
                    firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ) ).
                template cast< ObservationScalarType >( );

        // Compute spherical relative position for second transmitter wrt receiver.
        Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinatesSecondTransmitter =
                coordinate_conversions::convertCartesianToSpherical< ObservationScalarType >(
                    secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ) ).
                template cast< ObservationScalarType >( );

        // Compute right ascension and declination first transmitter.
        double rightAscensionFirstTransmitter = sphericalRelativeCoordinatesFirstTransmitter.z( );
        double declinationFirstTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesFirstTransmitter.y( );

        // Compute right ascension and declination second transmitter.
        double rightAscensionSecondTransmitter = sphericalRelativeCoordinatesSecondTransmitter.z( );
        double declinationSecondTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesSecondTransmitter.y( );

        // Compute partials of right ascension w.r.t. time.
        double timeDerivativeRightAscensionFirstTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                    firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                    firstTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );
        double timeDerivativeRightAscensionSecondTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                    secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                    secondTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );

        // Compute partials of declination w.r.t. time.
        double timeDerivativeDeclinationFirstTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                    firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                    firstTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );
        double timeDerivativeDeclinationSecondTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                    secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                    secondTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );

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
        double modifiedObservable = 1.0 / ( apparentDistance )
                * ( relativePositionInReceiverFrame[ 0 ] * relativeVelocityInReceiverFrame[ 0 ]
                + relativePositionInReceiverFrame[ 1 ] * relativeVelocityInReceiverFrame[ 1 ] );

        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << estimatedCentralInstant ).finished( );
    }


    //! Function to determine if the central instant observable is viable.
    bool isCentralInstantObservableViable( double time, double limitApparentDistance,
                                           std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations )
    {
        bool viableCentralInstant = false;
        if ( apparentDistanceObservations.size( ) > 3 )
        {
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > > apparentDistancesInterpolator
                    = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > >(
                        utilities::createVectorFromMapKeys< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ),
                        utilities::createVectorFromMapValues< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ), 4 );
            Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceAtObservationTime = apparentDistancesInterpolator->interpolate( time );
//        std::cout << "first apparent distance (is observable viable): " << apparentDistanceObservations.begin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "middle apparent distance (is observable viable): " << apparentDistanceAtObservationTime( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "final apparent distance (is observable viable): " << apparentDistanceObservations.rbegin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";

            if ( ( apparentDistanceObservations.begin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
                 && ( apparentDistanceObservations.rbegin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
                 && ( apparentDistanceAtObservationTime( 0, 0 ) <= limitApparentDistance ) )
            {
                viableCentralInstant = true;
            }
        }
        return viableCentralInstant;
    }


    //! Function to check whether the current central instant as already computed for a detected mutual approximation.
    bool isMutualApproximationAlreadyDetected( double estimatedCentralInstant )
    {
        bool mutualApproximationAlreadyDetected = false;
        for ( unsigned int i = 0 ; i < centralInstantsOfDetectedMutualApproximations_.size( ) ; i++ )
        {
            if ( ( estimatedCentralInstant - centralInstantsOfDetectedMutualApproximations_[ i ] ) / centralInstantsOfDetectedMutualApproximations_[ i ]
                 < 1.0e-6 )
            {
                mutualApproximationAlreadyDetected = true;
//                std::cout << "mutual approximation already detected!" << "\n\n";
            }
        }

        return mutualApproximationAlreadyDetected;
    }


    //! Function to get the object to calculate light time between the first transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the first transmitter and the receiver.
     * \return Object to calculate light time between the first transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return lightTimeCalculatorFirstTransmitter_;
    }

    //! Function to get the object to calculate light time between the second transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the second transmitter and the receiver.
     * \return Object to calculate light time between the second transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorSecondTransmitter( )
    {
        return lightTimeCalculatorSecondTransmitter_;
    }


    //! Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
    /*!
     * Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
     * \param List of viability calculators for apparent distance observables.
     */
    void setViabilityCalculatorsApparentDistances( std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances )
    {
        viabilityCalculatorsApparentDistances_.clear( );
        viabilityCalculatorsApparentDistances_ = viabilityCalculatorsApparentDistances;
    }


    //! Function to clear the former central instant observables already computed with this mutual approximation observation model.
    /*!
     * Function to clear the former central instant observables already computed with this mutual approximation observation model.
     */
    void clearAlreadyComputedCentralInstants( )
    {
        centralInstantsOfDetectedMutualApproximations_.clear( );
    }

    //! Function to retrieve the boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! while computing the central instant observable.
    bool isExistenceOfMutualApproximationChecked( )
    {
        return checkExistenceMutualApproximation_;
    }

    //! Function to retrieve the boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates.
    bool areObservableDuplicatesRevoved( )
    {
        return checkObservableDuplicates_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter_;

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter_;


    //! Apparent distance observation model object.
    std::shared_ptr< ApparentDistanceObservationModel< ObservationScalarType, TimeType > > apparentDistanceObservationModel_;

    //! Frequency at which the apparent distance between the two sources in the sky is collected.
    double frequencyApparentDistanceObservations_;

    //! Tolerance with respect to the expected central instant.
    double toleranceWrtExpectedCentralInstant_;

    //! Order of the polynomial fitting to derive the mutual approximation observable.
    int orderPolynomialFitting_;

    //! Upper limit value for the impact parameter, beyond which there is no close encounter.
    double upperLimitImpactParameter_;

    //! Boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! or not while computing the central instant observable (default value is true).
    bool checkExistenceMutualApproximation_;

    //! Boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates (default value is true).
    bool checkObservableDuplicates_;

    //! Settings objects for the root finder used to estimate the central instant t0;
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Vector storing the central instants of the already identified mutual approximations.
    std::vector< double > centralInstantsOfDetectedMutualApproximations_;

    //! List of observation viability calculators, which are used to discard simulated
    //! apparent distance if a given set of conditions are not fulfilled.
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances_;

};


//! Class for simulating mutual approximation observables, including impact paraneters.
/*!
 *  Class for simulating mutual approximation observables (including impact paraneters), using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class MutualApproximationWithImpactParameterObservationModel: public ObservationModel< 2, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     * Constructor
     * \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the first transmitter and the receiver.
     * \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the second transmitter and the receiver.
     * \param frequencyApparentDistanceObservations Time interval between two successive apparent distance measurements.
     * \param toleranceWrtExpectedCentralInstant Tolerance w.r.t. to the estimated central instant (time interval over which the central instant is looked for).
     * \param upperLimitImpactParameter Upper limit value for the impact parameter, beyond which there is no close encounter.
     * \param orderPolynomialFitting Order of the fitting polynomial to derive the central instant.
     * \param checkExistenceMutualApproximation Boolean denoting whether the existence of a mutual approximation at the given time must be verified.
     * \param checkObservableDuplicates Boolean denoting whether the duplicated central instant observables should be removed.
     * \param rootFinderSettings Root finder settings used to derive the central instant from the apparent distance observations.
     * \param mutualApproximationBiasCalculator Object for calculating system-dependent errors in the central instant
     *  observable, i.e. deviations from the physically ideal central instant observable between reference points (default none).
     * \param apparentDistancesBiasCalculator Object for calculating system-dependent errors in the apparent distance
     *  observables, i.e. deviations from the physically ideal apparent distance observables between reference points (default none).
     */
    MutualApproximationWithImpactParameterObservationModel(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter,
            const double frequencyApparentDistanceObservations,
            const double toleranceWrtExpectedCentralInstant,
            const double upperLimitImpactParameter,
            const int orderPolynomialFitting = 4,
            const bool checkExistenceMutualApproximation = true,
            const bool checkObservableDuplicates = true,
            const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
            const std::shared_ptr< ObservationBias< 2 > > mutualApproximationBiasCalculator = nullptr,
            const std::shared_ptr< ObservationBias< 1 > > apparentDistancesBiasCalculator = nullptr ):
        ObservationModel< 2, ObservationScalarType, TimeType >( mutual_approximation_with_impact_parameter, mutualApproximationBiasCalculator ),
        lightTimeCalculatorFirstTransmitter_( lightTimeCalculatorFirstTransmitter ),
        lightTimeCalculatorSecondTransmitter_( lightTimeCalculatorSecondTransmitter ),
        frequencyApparentDistanceObservations_( frequencyApparentDistanceObservations ),
        toleranceWrtExpectedCentralInstant_( toleranceWrtExpectedCentralInstant ),
        upperLimitImpactParameter_( upperLimitImpactParameter ),
        orderPolynomialFitting_( orderPolynomialFitting ),
        checkExistenceMutualApproximation_( checkExistenceMutualApproximation ),
        checkObservableDuplicates_( checkObservableDuplicates ),
        rootFinderSettings_( rootFinderSettings )
    {
        // Create required observation model to compute the intermediate apparent distance measurements.
        apparentDistanceObservationModel_ =
                std::make_shared< ApparentDistanceObservationModel< ObservationScalarType, TimeType > >
                ( lightTimeCalculatorFirstTransmitter_, lightTimeCalculatorSecondTransmitter_, apparentDistancesBiasCalculator );

        // Initialise empty vector of viability calculator for intermediate apparent distance observations.
        viabilityCalculatorsApparentDistances_ = std::vector< std::shared_ptr< ObservationViabilityCalculator > >( );

    }

    //! Destructor
    ~MutualApproximationWithImpactParameterObservationModel( ){ }

    //! Function to compute ideal central instants and impact parameters for mutual approximations.
    /*!
     *  This function compute ideal central instants and impact parameters for mutual approximations.
     *  The time argument should always be the reception time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which the mutual approximation is expected to occur.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value). Should always be receiver.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated central instant and impact parameter observable values.
     */
    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        std::cout.precision( 20 );

        double currentTime = time;

        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating mutual approximation observation, link end is not receiver." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState;


        ObservationScalarType estimatedCentralInstant = TUDAT_NAN;
        ObservationScalarType impactParameter = TUDAT_NAN;
        bool iterativeProcess = true;
        unsigned int iteration = 1;
        while ( iterativeProcess && ( iterativeProcess <= 2 ) ){

            if ( iteration == 2 )
            {
                iterativeProcess = false;
            }

        // Compute intermediate apparent distance observables.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations;

        double observationStartingTime = currentTime - toleranceWrtExpectedCentralInstant_;
        double observationEndingTime = currentTime + toleranceWrtExpectedCentralInstant_;

        for ( TimeType currentObservationTime = observationStartingTime ;
              currentObservationTime <= observationEndingTime ;
              currentObservationTime += frequencyApparentDistanceObservations_ )
        {
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceObservation =
                    apparentDistanceObservationModel_->computeObservationsWithLinkEndData( currentObservationTime, receiver, vectorOfTimes, vectorOfStates );

            bool observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, viabilityCalculatorsApparentDistances_ );
            if ( observationFeasible )
            {
                apparentDistanceObservations[ currentObservationTime ] = apparentDistanceObservation;
            }
        }


        /*ObservationScalarType*/ estimatedCentralInstant = TUDAT_NAN;
        /*ObservationScalarType*/ impactParameter = TUDAT_NAN;
        if ( ( !checkExistenceMutualApproximation_ ) || ( isCentralInstantObservableViable( currentTime, upperLimitImpactParameter_, apparentDistanceObservations ) ) )
        {
            // Polynomial fitting.
            std::vector< double > polynomialPowers;
            for ( unsigned int i = 0 ; i <= orderPolynomialFitting_ ; i++ )
            {
                polynomialPowers.push_back( i );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistances;
            apparentDistances.resize( apparentDistanceObservations.size( ) );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistanceTimes;
            apparentDistanceTimes.resize( apparentDistanceObservations.size( ) );

            double initialisationCounter = - ( int ) ( apparentDistanceObservations.size( ) / 2.0 );

            typename std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > >::iterator itr = apparentDistanceObservations.begin( );
            for ( unsigned int i = 0 ; i < apparentDistanceObservations.size( ) ; i++ )
            {
                apparentDistanceTimes( i, 0 ) = i + initialisationCounter; //itr->first / 3600.0;
                apparentDistances( i, 0 ) =  itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI;
                itr++;
            }

            // Compute coefficients of the fitting polynomial.
            Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceTimes, apparentDistances, polynomialPowers );

            Eigen::VectorXd derivativePolynomialCoefficients; derivativePolynomialCoefficients.resize( polynomialCoefficients.size( ) - 1 );
            for ( unsigned int i = 1 ; i <= orderPolynomialFitting_ ; i++ )
            {
                derivativePolynomialCoefficients[ i - 1 ] = i * polynomialCoefficients[ i ];
            }

            // Compute coefficients of the derivative of the fitting polynomial.
            double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
            double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
            double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];


            // Compute central instant analytically.
            double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
            double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

            double beta = Q * Q * Q + R * R;

            if ( beta < 0 )
            {
                double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
                double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
                double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
                double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "t2: " << t2 << "\n\n";
//                std::cout << "t3: " << t3 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t2 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t3 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
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
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
            }

//            std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";

            // Use of root finder given as input.
            std::shared_ptr< basic_mathematics::Function< double, double > > derivativePolynomialFittingFunction =
                    std::make_shared< PolynomialFunction >( derivativePolynomialCoefficients, orderPolynomialFitting_ - 1 );


            double estimatedCentralInstant2 = TUDAT_NAN;
            try
            {
                // Create root finder from root finder settings.
                std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder
                        = root_finders::createRootFinder( rootFinderSettings_, initialisationCounter, - initialisationCounter, 0.0 );
//                estimatedCentralInstant2 = observationStartingTime
//                    + ( rootFinder->execute( derivativePolynomialFittingFunction, 0.0 ) - initialisationCounter ) * frequencyApparentDistanceObservations_;

//                std::cout << "estimated central instant root finder: " << estimatedCentralInstant2 << "\n\n";
//                std::cout << "difference in central instant between analytical and root finder solutions: "
//                          << estimatedCentralInstant - estimatedCentralInstant2 << "\n\n";


                // Compute light-time and receiver/transmitter states.
                ObservationScalarType lightTimeFirstTransmitter = lightTimeCalculatorFirstTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, firstTransmitterState, estimatedCentralInstant, true );

                ObservationScalarType lightTimeSecondTransmitter = lightTimeCalculatorSecondTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, secondTransmitterState, estimatedCentralInstant, true );

                // Compute impact parameter at estimated central instant.
                impactParameter = apparentDistanceObservationModel_->computeIdealObservations( estimatedCentralInstant, receiver )[ 0 ];

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndStates.clear( );

                linkEndStates.push_back( firstTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( secondTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( receiverState.template cast< double >( ) );

                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeFirstTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeSecondTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant ) );

                if ( checkObservableDuplicates_ )
                {
                    if ( !isMutualApproximationAlreadyDetected( estimatedCentralInstant ) )
                    {
                        if ( !iterativeProcess )
                        {
                        centralInstantsOfDetectedMutualApproximations_.push_back( estimatedCentralInstant );
                        }
                    }
                    else
                    {
                        estimatedCentralInstant = TUDAT_NAN;
                        impactParameter = TUDAT_NAN;

                        // Set link end times and states.
                        linkEndTimes.clear( );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );

                        iterativeProcess = false;
                    }
                }
            }
            catch ( std::runtime_error& error )
            {
                std::cerr << "Warning, no mutual approximation is found around the requested observation time t = " << std::to_string( currentTime ) <<
                             ". The root finder has failed to find a minimum in the apparent distances history. Returns NAN as observable value." << std::endl;
                estimatedCentralInstant = TUDAT_NAN;
                impactParameter = TUDAT_NAN;

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );

                linkEndStates.clear( );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

                iterativeProcess = false;
            }

        }
        else
        {
            // Set link end times and states.
            linkEndTimes.clear( );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );

            linkEndStates.clear( );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

            iterativeProcess = false;
        }


        iteration++;
        currentTime = estimatedCentralInstant;

        }


        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << estimatedCentralInstant, impactParameter ).finished( );
    }


    //! Function to determine if the central instant observable is viable.
    bool isCentralInstantObservableViable( double time, double limitApparentDistance,
                                           std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations )
    {
        bool viableCentralInstant = false;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > > apparentDistancesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > >(
                    utilities::createVectorFromMapKeys< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ),
                    utilities::createVectorFromMapValues< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ), 4 );
        Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceAtObservationTime = apparentDistancesInterpolator->interpolate( time );
//        std::cout << "first apparent distance (is observable viable): " << apparentDistanceObservations.begin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "middle apparent distance (is observable viable): " << apparentDistanceAtObservationTime( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "final apparent distance (is observable viable): " << apparentDistanceObservations.rbegin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
        if ( ( apparentDistanceObservations.begin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceObservations.rbegin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceAtObservationTime( 0, 0 ) <= limitApparentDistance ) )
        {
            viableCentralInstant = true;
        }
        return viableCentralInstant;
    }


    //! Function to check whether the current central instant as already computed for a detected mutual approximation.
    bool isMutualApproximationAlreadyDetected( double estimatedCentralInstant )
    {
        bool mutualApproximationAlreadyDetected = false;
        for ( unsigned int i = 0 ; i < centralInstantsOfDetectedMutualApproximations_.size( ) ; i++ )
        {
            if ( ( estimatedCentralInstant - centralInstantsOfDetectedMutualApproximations_[ i ] ) / centralInstantsOfDetectedMutualApproximations_[ i ]
                 < 1.0e-6 )
            {
                mutualApproximationAlreadyDetected = true;
            }
        }

        return mutualApproximationAlreadyDetected;
    }


    //! Function to get the object to calculate light time between the first transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the first transmitter and the receiver.
     * \return Object to calculate light time between the first transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return lightTimeCalculatorFirstTransmitter_;
    }

    //! Function to get the object to calculate light time between the second transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the second transmitter and the receiver.
     * \return Object to calculate light time between the second transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorSecondTransmitter( )
    {
        return lightTimeCalculatorSecondTransmitter_;
    }


    //! Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
    /*!
     * Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
     * \param List of viability calculators for apparent distance observables.
     */
    void setViabilityCalculatorsApparentDistances( std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances )
    {
        viabilityCalculatorsApparentDistances_.clear( );
        viabilityCalculatorsApparentDistances_ = viabilityCalculatorsApparentDistances;
    }


    //! Function to clear the former central instant observables already computed with this mutual approximation observation model.
    /*!
     * Function to clear the former central instant observables already computed with this mutual approximation observation model.
     */
    void clearAlreadyComputedCentralInstants( )
    {
        centralInstantsOfDetectedMutualApproximations_.clear( );
    }

    //! Function to retrieve the boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! while computing the central instant observable.
    bool isExistenceOfMutualApproximationChecked( )
    {
        return checkExistenceMutualApproximation_;
    }

    //! Function to retrieve the boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates.
    bool areObservableDuplicatesRevoved( )
    {
        return checkObservableDuplicates_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter_;

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter_;


    //! Apparent distance observation model object.
    std::shared_ptr< ApparentDistanceObservationModel< ObservationScalarType, TimeType > > apparentDistanceObservationModel_;

    //! Frequency at which the apparent distance between the two sources in the sky is collected.
    double frequencyApparentDistanceObservations_;

    //! Tolerance with respect to the expected central instant.
    double toleranceWrtExpectedCentralInstant_;

    //! Order of the polynomial fitting to derive the mutual approximation observable.
    int orderPolynomialFitting_;

    //! Upper limit value for the impact parameter, beyond which there is no close encounter.
    double upperLimitImpactParameter_;

    //! Boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! or not while computing the central instant observable (default value is true).
    bool checkExistenceMutualApproximation_;

    //! Boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates (default value is true).
    bool checkObservableDuplicates_;

    //! Settings objects for the root finder used to estimate the central instant t0;
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Vector storing the central instants of the already identified mutual approximations.
    std::vector< double > centralInstantsOfDetectedMutualApproximations_;

    //! List of observation viability calculators, which are used to discard simulated
    //! apparent distance if a given set of conditions are not fulfilled.
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances_;

};



//! Class for simulating slightly modified mutual approximation observables.
/*!
 *  Class for simulating slightly modified mutual approximation observables (not central instants directly to simplify variational equations),
 *  using light-time (with light-time corrections) to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ModifiedMutualApproximationObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     * Constructor
     * \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the first transmitter and the receiver.
     * \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the second transmitter and the receiver.
     * \param frequencyApparentDistanceObservations Time interval between two successive apparent distance measurements.
     * \param toleranceWrtExpectedCentralInstant Tolerance w.r.t. to the estimated central instant (time interval over which the central instant is looked for).
     * \param upperLimitImpactParameter Upper limit value for the impact parameter, beyond which there is no close encounter.
     * \param orderPolynomialFitting Order of the fitting polynomial to derive the central instant.
     * \param checkExistenceMutualApproximation Boolean denoting whether the existence of a mutual approximation at the given time must be verified.
     * \param checkObservableDuplicates Boolean denoting whether the duplicated central instant observables should be removed.
     * \param rootFinderSettings Root finder settings used to derive the central instant from the apparent distance observations.
     * \param mutualApproximationBiasCalculator Object for calculating system-dependent errors in the central instant
     *  observable, i.e. deviations from the physically ideal central instant observable between reference points (default none).
     * \param apparentDistancesBiasCalculator Object for calculating system-dependent errors in the apparent distance
     *  observables, i.e. deviations from the physically ideal apparent distance observables between reference points (default none).
     */
    ModifiedMutualApproximationObservationModel(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter,
            const double frequencyApparentDistanceObservations,
            const double toleranceWrtExpectedCentralInstant,
            const double upperLimitImpactParameter,
            const int orderPolynomialFitting = 4,
            const bool checkExistenceMutualApproximation = true,
            const bool checkObservableDuplicates = true,
            const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
            const std::shared_ptr< ObservationBias< 1 > > mutualApproximationBiasCalculator = nullptr,
            const std::shared_ptr< ObservationBias< 1 > > apparentDistancesBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( mutual_approximation, mutualApproximationBiasCalculator ),
        lightTimeCalculatorFirstTransmitter_( lightTimeCalculatorFirstTransmitter ),
        lightTimeCalculatorSecondTransmitter_( lightTimeCalculatorSecondTransmitter ),
        frequencyApparentDistanceObservations_( frequencyApparentDistanceObservations ),
        toleranceWrtExpectedCentralInstant_( toleranceWrtExpectedCentralInstant ),
        upperLimitImpactParameter_( upperLimitImpactParameter ),
        orderPolynomialFitting_( orderPolynomialFitting ),
        checkExistenceMutualApproximation_( checkExistenceMutualApproximation ),
        checkObservableDuplicates_( checkObservableDuplicates ),
        rootFinderSettings_( rootFinderSettings )
    {
        // Create required observation model to compute the intermediate apparent distance measurements.
        apparentDistanceObservationModel_ =
                std::make_shared< ApparentDistanceObservationModel< ObservationScalarType, TimeType > >
                ( lightTimeCalculatorFirstTransmitter_, lightTimeCalculatorSecondTransmitter_, apparentDistancesBiasCalculator );

        // Initialise empty vector of viability calculator for intermediate apparent distance observations.
        viabilityCalculatorsApparentDistances_ = std::vector< std::shared_ptr< ObservationViabilityCalculator > >( );

    }

    //! Destructor
    ~ModifiedMutualApproximationObservationModel( ){ }

    //! Function to compute ideal slightly modified mutual approximation observables (not central instant directly to
    //! simplify the variational equations formulation).
    /*!
     *  Function to compute ideal slightly modified mutual approximation observables (not central instant directly to
     *  simplify the variational equations formulation).
     *  The time argument should always be the reception time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which the mutual approximation is expected to occur.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value). Should always be receiver.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated modified mutual approximation observable values.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        std::cout.precision( 20 );

        double currentTime = time;

        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating mutual approximation observation, link end is not receiver." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState;


        ObservationScalarType estimatedCentralInstant = TUDAT_NAN;
        ObservationScalarType modifiedMutualApproximationObservable = TUDAT_NAN;
        bool iterativeProcess = true;
        unsigned int iteration = 1;
        while ( iterativeProcess && ( iterativeProcess <= 2 ) ){

            if ( iteration == 2 )
            {
                iterativeProcess = false;
            }

        // Compute intermediate apparent distance observables.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations;

        double observationStartingTime = currentTime - toleranceWrtExpectedCentralInstant_;
        double observationEndingTime = currentTime + toleranceWrtExpectedCentralInstant_;

        // Compute apparent distances history.
        for ( TimeType currentObservationTime = observationStartingTime ;
              currentObservationTime <= observationEndingTime ;
              currentObservationTime += frequencyApparentDistanceObservations_ )
        {
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceObservation =
                    apparentDistanceObservationModel_->computeObservationsWithLinkEndData( currentObservationTime, receiver, vectorOfTimes, vectorOfStates );

            bool observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, viabilityCalculatorsApparentDistances_ );
            if ( observationFeasible )
            {
                apparentDistanceObservations[ currentObservationTime ] = apparentDistanceObservation;
            }
        }


        /*ObservationScalarType*/ estimatedCentralInstant = TUDAT_NAN;
        /*ObservationScalarType*/ modifiedMutualApproximationObservable = TUDAT_NAN;
        if ( ( !checkExistenceMutualApproximation_ ) || ( isCentralInstantObservableViable( currentTime, upperLimitImpactParameter_, apparentDistanceObservations ) ) )
        {

            // Polynomial fitting.
            std::vector< double > polynomialPowers;
            for ( unsigned int i = 0 ; i <= orderPolynomialFitting_ ; i++ )
            {
                polynomialPowers.push_back( i );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistances;
            apparentDistances.resize( apparentDistanceObservations.size( ) );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistanceTimes;
            apparentDistanceTimes.resize( apparentDistanceObservations.size( ) );

            double initialisationCounter = - ( int ) ( apparentDistanceObservations.size( ) / 2.0 );

            // Compute apparent
            typename std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > >::iterator itr = apparentDistanceObservations.begin( );
            for ( unsigned int i = 0 ; i < apparentDistanceObservations.size( ) ; i++ )
            {
                apparentDistanceTimes( i, 0 ) = i + initialisationCounter; //itr->first / 3600.0;
                apparentDistances( i, 0 ) =  itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI;
                itr++;
            }


            // Compute coefficients of the fitting polynomial.
            Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceTimes, apparentDistances, polynomialPowers );

            // Compute coefficients of the derivative of the fitting polynomial.
            Eigen::VectorXd derivativePolynomialCoefficients; derivativePolynomialCoefficients.resize( polynomialCoefficients.size( ) - 1 );
            for ( unsigned int i = 1 ; i <= orderPolynomialFitting_ ; i++ )
            {
                derivativePolynomialCoefficients[ i - 1 ] = i * polynomialCoefficients[ i ];
            }

            double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
            double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
            double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];



            // Compute central instant analytically.
            double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
            double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

            double beta = Q * Q * Q + R * R;

            if ( beta < 0 )
            {
                double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
                double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
                double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
                double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "t2: " << t2 << "\n\n";
//                std::cout << "t3: " << t3 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t2 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
                }
                if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t3 - initialisationCounter ) * frequencyApparentDistanceObservations_;
//                    std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";
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
//                std::cout << "t1: " << t1 << "\n\n";
//                std::cout << "initialisation counter: " << initialisationCounter << "\n\n";
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
            }

//            std::cout << "estimated central instant: " << estimatedCentralInstant << "\n\n";

            // Use of root finder given as input.
            std::shared_ptr< basic_mathematics::Function< double, double > > derivativePolynomialFittingFunction =
                    std::make_shared< PolynomialFunction >( derivativePolynomialCoefficients, orderPolynomialFitting_ - 1 );


            double estimatedCentralInstant2 = TUDAT_NAN;
            try
            {
                // Create root finder from root finder settings.
                std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder
                        = root_finders::createRootFinder( rootFinderSettings_, initialisationCounter, - initialisationCounter, 0.0 );
//                estimatedCentralInstant2 = observationStartingTime
//                    + ( rootFinder->execute( derivativePolynomialFittingFunction, 0.0 ) - initialisationCounter ) * frequencyApparentDistanceObservations_;

//                std::cout << "estimated central instant root finder: " << estimatedCentralInstant2 << "\n\n";
//                std::cout << "difference in central instant between analytical and root finder solutions: "
//                          << estimatedCentralInstant - estimatedCentralInstant2 << "\n\n";


                // Compute light-time and receiver/transmitter states.
                ObservationScalarType lightTimeFirstTransmitter = lightTimeCalculatorFirstTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, firstTransmitterState, estimatedCentralInstant, true );

                ObservationScalarType lightTimeSecondTransmitter = lightTimeCalculatorSecondTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, secondTransmitterState, estimatedCentralInstant, true );



                /// Compute slightly modified observable for different formulation of the variational equations.

                // Compute spherical relative position for first transmitter wrt receiver.
                Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinatesFirstTransmitter =
                        coordinate_conversions::convertCartesianToSpherical< ObservationScalarType >(
                            firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ) ).
                        template cast< ObservationScalarType >( );

                // Compute spherical relative position for second transmitter wrt receiver.
                Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinatesSecondTransmitter =
                        coordinate_conversions::convertCartesianToSpherical< ObservationScalarType >(
                            secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ) ).
                        template cast< ObservationScalarType >( );

                // Compute right ascension and declination first transmitter.
                double rightAscensionFirstTransmitter = sphericalRelativeCoordinatesFirstTransmitter.z( );
                double declinationFirstTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesFirstTransmitter.y( );

                // Compute right ascension and declination second transmitter.
                double rightAscensionSecondTransmitter = sphericalRelativeCoordinatesSecondTransmitter.z( );
                double declinationSecondTransmitter = mathematical_constants::PI / 2.0 - sphericalRelativeCoordinatesSecondTransmitter.y( );

                // Compute partials of right ascension w.r.t. time.
                double timeDerivativeRightAscensionFirstTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                            firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                            firstTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );
                double timeDerivativeRightAscensionSecondTransmitter = observation_partials::computePartialOfRightAscensionWrtTime(
                            secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                            secondTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );

                // Compute partials of declination w.r.t. time.
                double timeDerivativeDeclinationFirstTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                            firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                            firstTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );
                double timeDerivativeDeclinationSecondTransmitter = observation_partials::computePartialOfDeclinationWrtTime(
                            secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ),
                            secondTransmitterState.segment( 3, 3 ) - receiverState.segment( 3, 3 ) );

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
                modifiedMutualApproximationObservable = 1.0 / ( apparentDistance )
                        * ( relativePositionInReceiverFrame[ 0 ] * relativeVelocityInReceiverFrame[ 0 ]
                        + relativePositionInReceiverFrame[ 1 ] * relativeVelocityInReceiverFrame[ 1 ] );



                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndStates.clear( );

                linkEndStates.push_back( firstTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( secondTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( receiverState.template cast< double >( ) );

                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeFirstTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeSecondTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant ) );


                if ( checkObservableDuplicates_ )
                {
                    if ( !isMutualApproximationAlreadyDetected( estimatedCentralInstant ) )
                    {
                        if ( !iterativeProcess )
                        {
                        centralInstantsOfDetectedMutualApproximations_.push_back( estimatedCentralInstant );
                        }
                    }
                    else
                    {
                        estimatedCentralInstant = TUDAT_NAN;
                        modifiedMutualApproximationObservable = TUDAT_NAN;

                        // Set link end times and states.
                        linkEndTimes.clear( );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );

                        iterativeProcess = false;
                    }
                }
            }
            catch ( std::runtime_error& error )
            {
                std::cerr << "Warning, no mutual approximation is found around the requested observation time t = " << std::to_string( currentTime ) <<
                             ". The root finder has failed to find a minimum in the apparent distances history. Returns NAN as observable value." << std::endl;
                estimatedCentralInstant = TUDAT_NAN;
                modifiedMutualApproximationObservable = TUDAT_NAN;

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );

                linkEndStates.clear( );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

                iterativeProcess = false;
            }

        }
        else
        {
            // Set link end times and states.
            linkEndTimes.clear( );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );

            linkEndStates.clear( );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

            iterativeProcess = false;
        }


        iteration++;
        currentTime = estimatedCentralInstant;

        }


        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << modifiedMutualApproximationObservable ).finished( );
    }


    //! Function to determine if the central instant observable is viable.
    bool isCentralInstantObservableViable( double time, double limitApparentDistance,
                                           std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations )
    {
        bool viableCentralInstant = false;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > > apparentDistancesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > >(
                    utilities::createVectorFromMapKeys< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ),
                    utilities::createVectorFromMapValues< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ), 4 );
        Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceAtObservationTime = apparentDistancesInterpolator->interpolate( time );
//        std::cout << "first apparent distance (is observable viable): " << apparentDistanceObservations.begin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "middle apparent distance (is observable viable): " << apparentDistanceAtObservationTime( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
//        std::cout << "final apparent distance (is observable viable): " << apparentDistanceObservations.rbegin( )->second( 0, 0 ) * 3600.0 * 180.0 / mathematical_constants::PI << "\n\n";
        if ( ( apparentDistanceObservations.begin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceObservations.rbegin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceAtObservationTime( 0, 0 ) <= limitApparentDistance ) )
        {
            viableCentralInstant = true;
        }
        return viableCentralInstant;
    }


    //! Function to check whether the current central instant as already computed for a detected mutual approximation.
    bool isMutualApproximationAlreadyDetected( double estimatedCentralInstant )
    {
        bool mutualApproximationAlreadyDetected = false;
        for ( unsigned int i = 0 ; i < centralInstantsOfDetectedMutualApproximations_.size( ) ; i++ )
        {
            if ( ( estimatedCentralInstant - centralInstantsOfDetectedMutualApproximations_[ i ] ) / centralInstantsOfDetectedMutualApproximations_[ i ]
                 < 1.0e-6 )
            {
                mutualApproximationAlreadyDetected = true;
            }
        }

        return mutualApproximationAlreadyDetected;
    }


    //! Function to get the object to calculate light time between the first transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the first transmitter and the receiver.
     * \return Object to calculate light time between the first transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return lightTimeCalculatorFirstTransmitter_;
    }

    //! Function to get the object to calculate light time between the second transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the second transmitter and the receiver.
     * \return Object to calculate light time between the second transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorSecondTransmitter( )
    {
        return lightTimeCalculatorSecondTransmitter_;
    }


    //! Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
    /*!
     * Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
     * \param List of viability calculators for apparent distance observables.
     */
    void setViabilityCalculatorsApparentDistances( std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances )
    {
        viabilityCalculatorsApparentDistances_.clear( );
        viabilityCalculatorsApparentDistances_ = viabilityCalculatorsApparentDistances;
    }


    //! Function to clear the former central instant observables already computed with this mutual approximation observation model.
    /*!
     * Function to clear the former central instant observables already computed with this mutual approximation observation model.
     */
    void clearAlreadyComputedCentralInstants( )
    {
        centralInstantsOfDetectedMutualApproximations_.clear( );
    }

    //! Function to retrieve the boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! while computing the central instant observable.
    bool isExistenceOfMutualApproximationChecked( )
    {
        return checkExistenceMutualApproximation_;
    }

    //! Function to retrieve the boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates.
    bool areObservableDuplicatesRevoved( )
    {
        return checkObservableDuplicates_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter_;

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter_;


    //! Apparent distance observation model object.
    std::shared_ptr< ApparentDistanceObservationModel< ObservationScalarType, TimeType > > apparentDistanceObservationModel_;

    //! Frequency at which the apparent distance between the two sources in the sky is collected.
    double frequencyApparentDistanceObservations_;

    //! Tolerance with respect to the expected central instant.
    double toleranceWrtExpectedCentralInstant_;

    //! Order of the polynomial fitting to derive the mutual approximation observable.
    int orderPolynomialFitting_;

    //! Upper limit value for the impact parameter, beyond which there is no close encounter.
    double upperLimitImpactParameter_;

    //! Boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! or not while computing the central instant observable (default value is true).
    bool checkExistenceMutualApproximation_;

    //! Boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates (default value is true).
    bool checkObservableDuplicates_;

    //! Settings objects for the root finder used to estimate the central instant t0;
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Vector storing the central instants of the already identified mutual approximations.
    std::vector< double > centralInstantsOfDetectedMutualApproximations_;

    //! List of observation viability calculators, which are used to discard simulated
    //! apparent distance if a given set of conditions are not fulfilled.
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances_;

};



//! Class for simulating impact parameters as mutual approximation observables.
/*!
 *  Class for simulating impact parameters as mutual approximation observables, using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ImpactParameterMutualApproxObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     * Constructor
     * \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the first transmitter and the receiver.
     * \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the second transmitter and the receiver.
     * \param frequencyApparentDistanceObservations Time interval between two successive apparent distance measurements.
     * \param toleranceWrtExpectedCentralInstant Tolerance w.r.t. to the estimated central instant (time interval over which the central instant is looked for).
     * \param upperLimitImpactParameter Upper limit value for the impact parameter, beyond which there is no close encounter.
     * \param orderPolynomialFitting Order of the fitting polynomial to derive the central instant.
     * \param checkExistenceMutualApproximation Boolean denoting whether the existence of a mutual approximation at the given time must be verified.
     * \param checkObservableDuplicates Boolean denoting whether the duplicated central instant observables should be removed.
     * \param rootFinderSettings Root finder settings used to derive the central instant from the apparent distance observations.
     * \param mutualApproximationBiasCalculator Object for calculating system-dependent errors in the central instant
     *  observable, i.e. deviations from the physically ideal central instant observable between reference points (default none).
     * \param apparentDistancesBiasCalculator Object for calculating system-dependent errors in the apparent distance
     *  observables, i.e. deviations from the physically ideal apparent distance observables between reference points (default none).
     */
    ImpactParameterMutualApproxObservationModel(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter,
            const double frequencyApparentDistanceObservations,
            const double toleranceWrtExpectedCentralInstant,
            const double upperLimitImpactParameter,
            const int orderPolynomialFitting = 4,
            const bool checkExistenceMutualApproximation = true,
            const bool checkObservableDuplicates = true,
            const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-12, 60 ),
            const std::shared_ptr< ObservationBias< 1 > >impactParameterBiasCalculator = nullptr,
            const std::shared_ptr< ObservationBias< 1 > > apparentDistancesBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( impact_parameter_mutual_approx, impactParameterBiasCalculator ),
        lightTimeCalculatorFirstTransmitter_( lightTimeCalculatorFirstTransmitter ),
        lightTimeCalculatorSecondTransmitter_( lightTimeCalculatorSecondTransmitter ),
        frequencyApparentDistanceObservations_( frequencyApparentDistanceObservations ),
        toleranceWrtExpectedCentralInstant_( toleranceWrtExpectedCentralInstant ),
        upperLimitImpactParameter_( upperLimitImpactParameter ),
        orderPolynomialFitting_( orderPolynomialFitting ),
        checkExistenceMutualApproximation_( checkExistenceMutualApproximation ),
        checkObservableDuplicates_( checkObservableDuplicates ),
        rootFinderSettings_( rootFinderSettings )
    {
        // Create required observation model to compute the intermediate apparent distance measurements.
        apparentDistanceObservationModel_ =
                std::make_shared< ApparentDistanceObservationModel< ObservationScalarType, TimeType > >
                ( lightTimeCalculatorFirstTransmitter_, lightTimeCalculatorSecondTransmitter_, apparentDistancesBiasCalculator );

        // Initialise empty vector of viability calculator for intermediate apparent distance observations.
        viabilityCalculatorsApparentDistances_ = std::vector< std::shared_ptr< ObservationViabilityCalculator > >( );

    }

    //! Destructor
    ~ImpactParameterMutualApproxObservationModel( ){ }

    //! Function to compute ideal impact parameters for mutual approximations.
    /*!
     *  This function compute ideal impact parameters for mutual approximations.
     *  The time argument should always be the reception time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which the mutual approximation is expected to occur.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value). Should always be receiver.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated impact parameter observable values.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        std::cout.precision( 20 );

        double currentTime = time;

        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating impact parameter of mutual approximation, link end is not receiver." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState;


        ObservationScalarType estimatedCentralInstant = TUDAT_NAN;
        ObservationScalarType impactParameter = TUDAT_NAN;
        bool iterativeProcess = true;
        unsigned int iteration = 1;
        while ( iterativeProcess && ( iterativeProcess <= 2 ) ){

            if ( iteration == 2 )
            {
                iterativeProcess = false;
            }

        // Compute intermediate apparent distance observables.
        std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations;

        double observationStartingTime = currentTime - toleranceWrtExpectedCentralInstant_;
        double observationEndingTime = currentTime + toleranceWrtExpectedCentralInstant_;

        for ( TimeType currentObservationTime = observationStartingTime ;
              currentObservationTime <= observationEndingTime ;
              currentObservationTime += frequencyApparentDistanceObservations_ )
        {
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceObservation =
                    apparentDistanceObservationModel_->computeObservationsWithLinkEndData( currentObservationTime, receiver, vectorOfTimes, vectorOfStates );

            bool observationFeasible = isObservationViable( vectorOfStates, vectorOfTimes, viabilityCalculatorsApparentDistances_ );
            if ( observationFeasible )
            {
                apparentDistanceObservations[ currentObservationTime ] = apparentDistanceObservation;
            }
        }


        /*ObservationScalarType*/ estimatedCentralInstant = TUDAT_NAN;
        /*ObservationScalarType*/ impactParameter = TUDAT_NAN;
        if ( ( !checkExistenceMutualApproximation_ ) || ( isCentralInstantObservableViable( currentTime, upperLimitImpactParameter_, apparentDistanceObservations ) ) )
        {
            // Polynomial fitting.
            std::vector< double > polynomialPowers;
            for ( unsigned int i = 0 ; i <= orderPolynomialFitting_ ; i++ )
            {
                polynomialPowers.push_back( i );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistances;
            apparentDistances.resize( apparentDistanceObservations.size( ) );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > apparentDistanceTimes;
            apparentDistanceTimes.resize( apparentDistanceObservations.size( ) );

            double initialisationCounter = - ( int ) ( apparentDistanceObservations.size( ) / 2.0 );

            typename std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > >::iterator itr = apparentDistanceObservations.begin( );
            for ( unsigned int i = 0 ; i < apparentDistanceObservations.size( ) ; i++ )
            {
                apparentDistanceTimes( i, 0 ) = i + initialisationCounter;
                apparentDistances( i, 0 ) =  itr->second[ 0 ] * 3600.0 * 180.0 / mathematical_constants::PI;
                itr++;
            }

            // Compute coefficients of the fitting polynomial.
            Eigen::VectorXd polynomialCoefficients = linear_algebra::getLeastSquaresPolynomialFit( apparentDistanceTimes, apparentDistances, polynomialPowers );

            Eigen::VectorXd derivativePolynomialCoefficients; derivativePolynomialCoefficients.resize( polynomialCoefficients.size( ) - 1 );
            for ( unsigned int i = 1 ; i <= orderPolynomialFitting_ ; i++ )
            {
                derivativePolynomialCoefficients[ i - 1 ] = i * polynomialCoefficients[ i ];
            }

            // Compute coefficients of the derivative of the fitting polynomial.
            double b0 = derivativePolynomialCoefficients[ 0 ] / derivativePolynomialCoefficients[ 3 ];
            double b1 = derivativePolynomialCoefficients[ 1 ] / derivativePolynomialCoefficients[ 3 ];
            double b2 = derivativePolynomialCoefficients[ 2 ] / derivativePolynomialCoefficients[ 3 ];


            // Compute central instant analytically.
            double Q = ( 3.0 * b1 - b2 * b2 ) / 9.0;
            double R = ( 9.0 * b2 * b1 - 27.0 * b0 - 2.0 * b2 * b2 * b2 ) / 54.0;

            double beta = Q * Q * Q + R * R;

            if ( beta < 0 )
            {
                double theta = std::acos( R / std::sqrt( - Q * Q * Q ) );
                double t1 = 2.0 * std::sqrt( - Q ) * std::cos( theta / 3.0 ) - b2 / 3.0;
                double t2 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta  + 2.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
                double t3 = 2.0 * std::sqrt( - Q ) * std::cos( ( theta + 4.0 * mathematical_constants::PI ) / 3.0 ) - b2 / 3.0;
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
                if ( t2 >= initialisationCounter && t2 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t2 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
                if ( t3 >= initialisationCounter && t3 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t3 - initialisationCounter ) * frequencyApparentDistanceObservations_;
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
                if ( t1 >= initialisationCounter && t1 <= - initialisationCounter )
                {
                    estimatedCentralInstant = observationStartingTime + ( t1 - initialisationCounter ) * frequencyApparentDistanceObservations_;
                }
            }

            // Use of root finder given as input.
            std::shared_ptr< basic_mathematics::Function< double, double > > derivativePolynomialFittingFunction =
                    std::make_shared< PolynomialFunction >( derivativePolynomialCoefficients, orderPolynomialFitting_ - 1 );


            double estimatedCentralInstant2 = TUDAT_NAN;
            try
            {
                // Create root finder from root finder settings.
                std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder
                        = root_finders::createRootFinder( rootFinderSettings_, initialisationCounter, - initialisationCounter, 0.0 );
//                estimatedCentralInstant2 = observationStartingTime
//                    + ( rootFinder->execute( derivativePolynomialFittingFunction, 0.0 ) - initialisationCounter ) * frequencyApparentDistanceObservations_;

                // Compute light-time and receiver/transmitter states.
                ObservationScalarType lightTimeFirstTransmitter = lightTimeCalculatorFirstTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, firstTransmitterState, estimatedCentralInstant, true );

                ObservationScalarType lightTimeSecondTransmitter = lightTimeCalculatorSecondTransmitter_->calculateLightTimeWithLinkEndsStates(
                            receiverState, secondTransmitterState, estimatedCentralInstant, true );

                // Compute impact parameter at estimated central instant.
                impactParameter = apparentDistanceObservationModel_->computeIdealObservations( estimatedCentralInstant, receiver )[ 0 ];

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndStates.clear( );

                linkEndStates.push_back( firstTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( secondTransmitterState.template cast< double >( ) );
                linkEndStates.push_back( receiverState.template cast< double >( ) );

                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeFirstTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant - lightTimeSecondTransmitter ) );
                linkEndTimes.push_back( static_cast< double >( estimatedCentralInstant ) );

                if ( checkObservableDuplicates_ )
                {
                    if ( !isMutualApproximationAlreadyDetected( estimatedCentralInstant ) )
                    {
                        if ( !iterativeProcess )
                        {
                        centralInstantsOfDetectedMutualApproximations_.push_back( estimatedCentralInstant );
                        }
                    }
                    else
                    {
                        estimatedCentralInstant = TUDAT_NAN;
                        impactParameter = TUDAT_NAN;

                        // Set link end times and states.
                        linkEndTimes.clear( );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );
                        linkEndTimes.push_back( TUDAT_NAN );

                        iterativeProcess = false;
                    }
                }
            }
            catch ( std::runtime_error& error )
            {
                std::cerr << "Warning, no mutual approximation is found around the requested observation time t = " << std::to_string( currentTime ) <<
                             ". The root finder has failed to find a minimum in the apparent distances history. Returns NAN as observable value." << std::endl;
                estimatedCentralInstant = TUDAT_NAN;
                impactParameter = TUDAT_NAN;

                // Set link end times and states.
                linkEndTimes.clear( );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );
                linkEndTimes.push_back( TUDAT_NAN );

                linkEndStates.clear( );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
                linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

                iterativeProcess = false;
            }

        }
        else
        {
            // Set link end times and states.
            linkEndTimes.clear( );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );
            linkEndTimes.push_back( TUDAT_NAN );

            linkEndStates.clear( );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );
            linkEndStates.push_back( TUDAT_NAN * Eigen::Vector6d::Ones( ) );

            iterativeProcess = false;
        }


        iteration++;
        currentTime = estimatedCentralInstant;

        }


        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << impactParameter ).finished( );
    }


    //! Function to determine if the central instant observable is viable.
    bool isCentralInstantObservableViable( double time, double limitApparentDistance,
                                           std::map< TimeType, Eigen::Matrix< ObservationScalarType, 1, 1 > > apparentDistanceObservations )
    {
        bool viableCentralInstant = false;
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > > apparentDistancesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< ObservationScalarType, 1, 1 > > >(
                    utilities::createVectorFromMapKeys< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ),
                    utilities::createVectorFromMapValues< Eigen::Matrix< ObservationScalarType, 1, 1 >, double >( apparentDistanceObservations ), 4 );
        Eigen::Matrix< ObservationScalarType, 1, 1 > apparentDistanceAtObservationTime = apparentDistancesInterpolator->interpolate( time );
        if ( ( apparentDistanceObservations.begin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceObservations.rbegin( )->second( 0, 0 ) >= apparentDistanceAtObservationTime( 0, 0 ) )
             && ( apparentDistanceAtObservationTime( 0, 0 ) <= limitApparentDistance ) )
        {
            viableCentralInstant = true;
        }
        return viableCentralInstant;
    }


    //! Function to check whether the current central instant as already computed for a detected mutual approximation.
    bool isMutualApproximationAlreadyDetected( double estimatedCentralInstant )
    {
        bool mutualApproximationAlreadyDetected = false;
        for ( unsigned int i = 0 ; i < centralInstantsOfDetectedMutualApproximations_.size( ) ; i++ )
        {
            if ( ( estimatedCentralInstant - centralInstantsOfDetectedMutualApproximations_[ i ] ) / centralInstantsOfDetectedMutualApproximations_[ i ]
                 < 1.0e-6 )
            {
                mutualApproximationAlreadyDetected = true;
            }
        }

        return mutualApproximationAlreadyDetected;
    }


    //! Function to get the object to calculate light time between the first transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the first transmitter and the receiver.
     * \return Object to calculate light time between the first transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return lightTimeCalculatorFirstTransmitter_;
    }

    //! Function to get the object to calculate light time between the second transmitter and the receiver.
    /*!
     * Function to get the object to calculate light time between the second transmitter and the receiver.
     * \return Object to calculate light time between the second transmitter and the receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorSecondTransmitter( )
    {
        return lightTimeCalculatorSecondTransmitter_;
    }


    //! Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
    /*!
     * Function to set the viability calculators for the apparent distance measurements used to derive the central instant
    //! of the mutual approximation.
     * \param List of viability calculators for apparent distance observables.
     */
    void setViabilityCalculatorsApparentDistances( std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances )
    {
        viabilityCalculatorsApparentDistances_.clear( );
        viabilityCalculatorsApparentDistances_ = viabilityCalculatorsApparentDistances;
    }


    //! Function to clear the former central instant observables already computed with this mutual approximation observation model.
    /*!
     * Function to clear the former central instant observables already computed with this mutual approximation observation model.
     */
    void clearAlreadyComputedCentralInstants( )
    {
        centralInstantsOfDetectedMutualApproximations_.clear( );
    }

    //! Function to retrieve the boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! while computing the central instant observable.
    bool isExistenceOfMutualApproximationChecked( )
    {
        return checkExistenceMutualApproximation_;
    }

    //! Function to retrieve the boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates.
    bool areObservableDuplicatesRevoved( )
    {
        return checkObservableDuplicates_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter_;

    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter_;


    //! Apparent distance observation model object.
    std::shared_ptr< ApparentDistanceObservationModel< ObservationScalarType, TimeType > > apparentDistanceObservationModel_;

    //! Frequency at which the apparent distance between the two sources in the sky is collected.
    double frequencyApparentDistanceObservations_;

    //! Tolerance with respect to the expected central instant.
    double toleranceWrtExpectedCentralInstant_;

    //! Order of the polynomial fitting to derive the mutual approximation observable.
    int orderPolynomialFitting_;

    //! Upper limit value for the impact parameter, beyond which there is no close encounter.
    double upperLimitImpactParameter_;

    //! Boolean denoting whether the existence of a mutual approximation at the given time must be verified
    //! or not while computing the central instant observable (default value is true).
    bool checkExistenceMutualApproximation_;

    //! Boolean denoting whether the current central instant observable should be compared with former observables already computed with
    //! the mutual approximation observation model, in order to remove duplicates (default value is true).
    bool checkObservableDuplicates_;

    //! Settings objects for the root finder used to estimate the central instant t0;
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Vector storing the central instants of the already identified mutual approximations.
    std::vector< double > centralInstantsOfDetectedMutualApproximations_;

    //! List of observation viability calculators, which are used to discard simulated
    //! apparent distance if a given set of conditions are not fulfilled.
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculatorsApparentDistances_;

};



} // namespace observation_models

} // namespace tudat

#endif // TUDAT_MUTUALAPPROXIMATIONOBSERVATIONMODEL_H
