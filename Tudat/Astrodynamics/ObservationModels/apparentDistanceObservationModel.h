/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_APPARENTDISTANCEOBSERVATIONMODEL_H
#define TUDAT_APPARENTDISTANCEOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating apparent distance observables (between two objects in the sky).
/*!
 *  Class for simulating apparent distance observables (between two objects in the sky), using light-time (with light-time corrections)
 *  to determine the states of the link ends (sources and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class ApparentDistanceObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the first transmitter and the receiver.
     *  \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between the second transmitter and the receiver.
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    ApparentDistanceObservationModel(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( apparent_distance, observationBiasCalculator ),
        lightTimeCalculatorFirstTransmitter_( lightTimeCalculatorFirstTransmitter ),
        lightTimeCalculatorSecondTransmitter_( lightTimeCalculatorSecondTransmitter ){ }

    //! Destructor
    ~ApparentDistanceObservationModel( ){ }

    //! Function to compute ideal apparent position between two objects at given time.
    /*!
     *  This function compute ideal apparent position between two objects at given time. The time argument should always be the
     *  reception time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value). It should be the reception time, otherwise an error is returned.
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated apparent distance observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating apparent distance observation, link end is not receiver." );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState;

        // Compute light-time and receiver/transmitter states.
        ObservationScalarType lightTimeFirstTransmitter = lightTimeCalculatorFirstTransmitter_->calculateLightTimeWithLinkEndsStates(
                    receiverState, firstTransmitterState, time, true );

        ObservationScalarType lightTimeSecondTransmitter = lightTimeCalculatorSecondTransmitter_->calculateLightTimeWithLinkEndsStates(
                    receiverState, secondTransmitterState, time, true );

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

        // Compute average declination.
        double averageDeclination = ( declinationFirstTransmitter + declinationSecondTransmitter ) / 2.0;

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndStates.push_back( firstTransmitterState.template cast< double >( ) );
        linkEndStates.push_back( secondTransmitterState.template cast< double >( ) );
        linkEndStates.push_back( receiverState.template cast< double >( ) );

        linkEndTimes.push_back( static_cast< double >( time - lightTimeFirstTransmitter ) );
        linkEndTimes.push_back( static_cast< double >( time - lightTimeSecondTransmitter ) );
        linkEndTimes.push_back( static_cast< double >( time ) );

        // Return observable.
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) <<
                 std::sqrt( ( ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter )
                              * std::cos( averageDeclination ) )
                            * ( ( rightAscensionSecondTransmitter - rightAscensionFirstTransmitter )
                                * std::cos( averageDeclination ) )
                            + ( declinationSecondTransmitter - declinationFirstTransmitter )
                            * ( declinationSecondTransmitter - declinationFirstTransmitter ) ) ).finished( );
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


private:

    //! Object to calculate light time between the first transmitter and the receiver.
    /*!
     *  Object to calculate light time between the first transmitter and the receiver, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter_;

    //! Object to calculate light time between the second transmitter and the receiver.
    /*!
     *  Object to calculate light time between the second transmitter and the receiver, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_APPARENTDISTANCEOBSERVATIONMODEL_H
