/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/apparentDistancePartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the right ascension and declination of the observed body as seen from the observer.
std::pair< double, double > computeRightAscensionAndDeclination( Eigen::Vector3d relativePosition )
{
    Eigen::Matrix< double, 3, 1 > sphericalCoordinates = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                - relativePosition ).template cast< double >( );
    double rightAscension = sphericalCoordinates.z( );
    double declination = mathematical_constants::PI / 2.0 - sphericalCoordinates.y( );

    return std::make_pair( rightAscension, declination );
}


//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - relativeRangeVector( 1 );
    partial( 1 ) = relativeRangeVector( 0 );
    partial /= ( partialMultiplier * ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) +
                                       relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    double range = relativeRangeVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - relativeRangeVector( 0 ) * relativeRangeVector( 2 );
    partial( 1 ) = - relativeRangeVector( 1 ) * relativeRangeVector( 2 );
    partial( 2 ) = relativeRangeVector( 0 ) * relativeRangeVector( 0 ) + relativeRangeVector( 1 ) * relativeRangeVector( 1 );
    partial *= partialMultiplier /
            ( range * range * std::sqrt( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) + relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );

    return partial;
}

////! Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
////! observed object.
//Eigen::Matrix< double, 2, 3 > computePartialOfAngularPositionWrtLinkEndPosition(
//        Eigen::Vector3d relativeRangeVector,
//        const bool isLinkEndReceiver )
//{
//    Eigen::Matrix< double, 2, 3 > angularPositionPartial;
//    angularPositionPartial.block( 0, 0, 1, 3 ) = computePartialOfRightAscensionWrtLinkEndPosition(
//                relativeRangeVector, isLinkEndReceiver );
//    angularPositionPartial.block( 1, 0, 1, 3 ) = computePartialOfDeclinationWrtLinkEndPosition(
//                relativeRangeVector, isLinkEndReceiver );
//    return angularPositionPartial;
//}

//! Update the scaling object to the current times and states
void ApparentDistanceScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                      const std::vector< double >& times,
                                      const observation_models::LinkEndType fixedLinkEnd,
                                      const observation_models::LinkEnds linkEnds,
                                      const Eigen::VectorXd currentObservation )
{
    if ( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating the apparent distance scaling object, fixed link end time different from receiver." );
    }

    Eigen::Vector6d firstTransmitterState = linkEndStates[ 0 ];
    Eigen::Vector6d secondTransmitterState = linkEndStates[ 1 ];
    Eigen::Vector6d receiverState = linkEndStates[ 2 ];

    Eigen::Vector3d relativeRangeVectorFirstTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRangeVectorSecondTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

//    Eigen::Vector3d normalizedRelativeRangeVectorFirstTransmitter = relativeRangeVectorFirstTransmitter.normalized( );
//    Eigen::Vector3d normalizedRelativeRangeVectorSecondTransmitter = relativeRangeVectorSecondTransmitter.normalized( );

    std::pair< double, double > rightAscensionAndDeclinationFirstTransmitter =
            computeRightAscensionAndDeclination( relativeRangeVectorFirstTransmitter );
    std::pair< double, double > rightAscensionAndDeclinationSecondTransmitter =
            computeRightAscensionAndDeclination( relativeRangeVectorSecondTransmitter );

    rightAscensionFirstTransmitter_ = rightAscensionAndDeclinationFirstTransmitter.first;
    declinationFirstTransmitter_ = rightAscensionAndDeclinationFirstTransmitter.second;
    rightAscensionSecondTransmitter_ = rightAscensionAndDeclinationSecondTransmitter.first;
    declinationSecondTransmitter_ = rightAscensionAndDeclinationSecondTransmitter.second;

    averageDeclination_ = ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0;

    XcoordinateReceiverFrame_ = ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
            * std::cos( averageDeclination_ );
    YcoordinateReceiverFrame_ = declinationSecondTransmitter_ - declinationFirstTransmitter_;


    // Compute reference scaling factors for angular position of both transmitters (at receiver).
    std::vector< Eigen::Vector6d > linkEndStatesFirstTransmitterReceiver;
    linkEndStatesFirstTransmitterReceiver.push_back( firstTransmitterState );
    linkEndStatesFirstTransmitterReceiver.push_back( receiverState );

    std::vector< Eigen::Vector6d > linkEndStatesSecondTransmitterReceiver;
    linkEndStatesSecondTransmitterReceiver.push_back( secondTransmitterState );
    linkEndStatesSecondTransmitterReceiver.push_back( receiverState );

    referenceAngularPositionScaler_->update( linkEndStatesFirstTransmitterReceiver, times, fixedLinkEnd, linkEnds, currentObservation );
    angularPositionScalingFactorFirstTransmitter_ = referenceAngularPositionScaler_->getScalingFactor( observation_models::receiver );
    angularPositionLightTimeCorrectionScalingFirstTransmitter_ = referenceAngularPositionScaler_->getLightTimePartialScalingFactor( );

    referenceAngularPositionScaler_->update( linkEndStatesSecondTransmitterReceiver, times, fixedLinkEnd, linkEnds, currentObservation );
    angularPositionScalingFactorSecondTransmitter_ = referenceAngularPositionScaler_->getScalingFactor( observation_models::receiver );
    angularPositionLightTimeCorrectionScalingSecondTransmitter_ = referenceAngularPositionScaler_->getLightTimePartialScalingFactor( );


    //


//    // Compute common scaling factor
//    scalingFactor_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVectorFirstTransmitter, true );




//    // Compute scaling for receiver reference
//    if( fixedLinkEnd == observation_models::receiver )
//    {
//        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
//                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
//        referenceScalingFactor_ =
//                scalingFactor_ *
//                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
//                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
//    }
//    // Compute scaling for transmitter reference
//    else if( fixedLinkEnd == observation_models::transmitter )
//    {
//        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
//                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
//        referenceScalingFactor_ =
//                scalingFactor_ *
//                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
//                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
//    }


    // TO BE MODIFIED, WORKS AS WELL FOR PARTIALS WRT LIGHT TIME CORRECTION PARAMETERS
    scalingFactorXCoefficient_ = std::make_pair(
                std::cos( averageDeclination_ ), - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
                * std::sin( averageDeclination_ ) );
    scalingFactorYCoefficient_ = 1.0;

    XscalingFactors_ = //scalingFactorsXrightAscensionContribution_ =
            std::make_pair( - angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) * scalingFactorXCoefficient_.first
                            + angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorXCoefficient_.second,
                            angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 ) * scalingFactorXCoefficient_.first
                            + angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorXCoefficient_.second );

    XlightTimeCorrectionScalingFactors_ = //scalingFactorsXrightAscensionContribution_ =
            std::make_pair( - angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 0 ] * scalingFactorXCoefficient_.first
                            + angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 1 ] * scalingFactorXCoefficient_.second,
                            angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 0 ] * scalingFactorXCoefficient_.first
                            + angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 1 ] * scalingFactorXCoefficient_.second );
//    partialsXrightAscensionContribution =
//            ( angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 ) - angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) )
//            * std::cos( averageDeclination_ );

//    scalingFactorsXdeclinationContribution_ =
//            std::make_pair( - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
//                            * std::sin( averageDeclination_ ) * angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ),
//                            - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
//                            * std::sin( averageDeclination_ ) * angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) );
////    partialsXdeclinationContribution *=  - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
////            * std::sin( averageDeclination_ );
////            - ( rightAscensionAndDeclinationSecondTransmitter.first - rightAscensionAndDeclinationFirstTransmitter.first ) / 2.0
////            * std::sin( averageDeclination_ )
////            * ( angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) + angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) );

    YscalingFactors_ =
            std::make_pair( - angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorYCoefficient_,
            angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorYCoefficient_ );

    YlightTimeCorrectionScalingFactors_ =
            std::make_pair( - angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 1 ] * scalingFactorYCoefficient_,
            angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 1 ] * scalingFactorYCoefficient_ );

    //Eigen::Matrix< double, 1, 3 >::Zero( 1, 3 );
//    partialsY = angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) - angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 );

//    double X = std::cos( averageDeclination_ )
//            * ( rightAscensionAndDeclinationSecondTransmitter.first - rightAscensionAndDeclinationFirstTransmitter.first );

//    double Y = rightAscensionAndDeclinationSecondTransmitter.second - rightAscensionAndDeclinationFirstTransmitter.second;

//    referenceScalingFactor_ = 1.0 / std::sqrt( X * X + Y * Y )
//            * ( X * partialsXWrtLinkEndPosition + Y * partialsY );

//    // UPDATE REFERENCE LIGHT TIME CORRECTION SCALING

    currentLinkEndType_ = fixedLinkEnd;

}

//! Function to calculate the observation partial(s) at required time and state
ApparentDistancePartial::ApparentDistancePartialReturnType ApparentDistancePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    if( linkEndOfFixedTime != apparentDistanceScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error apparent distance partial and scaling are inconsistent" );
    }

    ApparentDistancePartialReturnType returnPartial;

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        if( positionPartialIterator_->first == observation_models::transmitter )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];
        }
        else if( positionPartialIterator_->first == observation_models::transmitter2 )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];
        }
        else if( positionPartialIterator_->first == observation_models::receiver )
        {
            currentState_ = states[ 2 ];
            currentTime_ = times[ 2 ];
        }

        // Scale position partials
        returnPartial.push_back(
                    std::make_pair(
                        apparentDistanceScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }


    // Add scaled light-time correcion partials.
    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ); i++ )
    {
        currentLinkTimeCorrectionPartialFirstTransmitter_ = lighTimeCorrectionPartialsFunctionsFirstTransmitter_.at( i )( states, times );
        currentLinkTimeCorrectionPartialSecondTransmitter_ = lighTimeCorrectionPartialsFunctionsSecondTransmitter_.at( i )( states, times );

        if ( currentLinkTimeCorrectionPartialFirstTransmitter_.second !=
             currentLinkTimeCorrectionPartialSecondTransmitter_.second )
        {
            throw std::runtime_error( "Error when making apparent distance light time correction partials, unconsistency"
                                      " in receiver times between receiver - first transmitter and receiver - second transmitter legs." );
        }

        returnPartial.push_back(
                    std::make_pair( apparentDistanceScaler_->getLightTimePartialScalingFactor( ) *
                                    physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartialFirstTransmitter_.first,
                    currentLinkTimeCorrectionPartialFirstTransmitter_.second ) );
    }


    return returnPartial;
}

}

}
