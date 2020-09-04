/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/ObservationModels/observationViabilityCalculator.h"

namespace tudat
{

namespace observation_models
{

//! Function to check whether an observation is viable
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times, const LinkEnds& linkEnds,
        const std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > >& viabilityCalculators )
{
    bool isObservationFeasible = 1;

    if( viabilityCalculators.count( linkEnds ) > 0 )
    {
        isObservationFeasible = isObservationViable( states, times, viabilityCalculators.at( linkEnds ) );
    }

    return isObservationFeasible;
}

//! Function to check whether an observation is viable
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times,
        const std::vector< std::shared_ptr< ObservationViabilityCalculator > >& viabilityCalculators )
{
    bool isObservationFeasible = 1;

    for( unsigned int i = 0; i < viabilityCalculators.size( ); i++ )
    {
        if( viabilityCalculators.at( i )->isObservationViable( states, times ) == 0 )
        {
            isObservationFeasible = 0;
            break;
        }
    }

    return isObservationFeasible;
}

//! Function for determining whether the elevation angle at station is sufficient to allow observation
bool MinimumElevationAngleCalculator::isObservationViable(
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;

    // Iterate over all sets of entries of input vector for which elvation angle is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        // Check if elevation angle criteria is met for current link.
        if( ground_stations::isTargetInView(
                    linkEndTimes.at( linkEndIndices_.at( i ).first ),
                    ( linkEndStates.at( linkEndIndices_.at( i ).second ) - linkEndStates.at( linkEndIndices_.at( i ).first ) )
                    .segment( 0, 3 ), pointingAngleCalculator_, minimumElevationAngle_ ) == 0 )
        {
            isObservationPossible = 0;
        }
    }

    return isObservationPossible;
}

//! Function for determining whether the avoidance angle to a given body at station is sufficient to allow observation.
bool BodyAvoidanceAngleCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                        const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;
    Eigen::Vector3d positionOfBodyToAvoid;
    double currentCosineOfAngle;

    // Iterate over all sets of entries of input vector for which avoidance angle is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        // Compute cosine of avoidance angles
        positionOfBodyToAvoid = stateFunctionOfBodyToAvoid_(
                    ( linkEndTimes.at( linkEndIndices_.at( i ).first ) + linkEndTimes.at( linkEndIndices_.at( i ).second ) ) / 2.0 )
                .segment( 0, 3 );
        currentCosineOfAngle = linear_algebra::computeCosineOfAngleBetweenVectors(
                    positionOfBodyToAvoid - ( linkEndStates.at( linkEndIndices_.at( i ).first ) ).segment( 0, 3 ),
                    linkEndStates.at( linkEndIndices_.at( i ).second ).segment( 0, 3 ) -
                    linkEndStates.at( linkEndIndices_.at( i ).first ).segment( 0, 3 ) );

        // Check if avoidance angle is sufficiently large
        if( currentCosineOfAngle > std::cos( bodyAvoidanceAngle_ ) )
        {
            isObservationPossible = 0;
            break;
        }
    }

    return isObservationPossible;
}

//! Function for determining whether the link is occulted during the observataion.
bool OccultationCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                 const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;
    Eigen::Vector3d positionOfOccultingBody;


    // Iterate over all sets of entries of input vector for which occultation is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        // Get position of occulting body
        positionOfOccultingBody = stateFunctionOfOccultingBody_(
                    ( linkEndTimes.at( linkEndIndices_.at( i ).first ) +
                      linkEndTimes.at( linkEndIndices_.at( i ).second ) ) / 2.0 ).segment( 0, 3 );

        // Check if observing link end is occulted by body.
        if( mission_geometry::computeShadowFunction(
                    linkEndStates.at( linkEndIndices_.at( i ).first ).segment( 0, 3 ), 0.0,
                    positionOfOccultingBody,
                    radiusOfOccultingBody_,
                    linkEndStates.at( linkEndIndices_.at( i ).second ).segment( 0, 3 ) ) < 1.0E-10 )
        {
            isObservationPossible = 0;
            break;
        }
    }//

    return isObservationPossible;
}


//! Function for determining whether there is a mutual approximation and thus whether the central instant is defined or not.
bool MutualApproximationCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                         const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = false;
//    Eigen::Vector3d positionOfOccultingBody;


    // Iterate over all sets of entries of input vector for which occultation is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        if ( !isnan( linkEndTimes.at( linkEndIndices_.at( i ).first ) )
             || !isnan( linkEndTimes.at( linkEndIndices_.at( i ).second ) ) )
        {
//        // Get position of occulting body
//        positionOfOccultingBody = stateFunctionOfOccultingBody_(
//                    ( linkEndTimes.at( linkEndIndices_.at( i ).first ) +
//                      linkEndTimes.at( linkEndIndices_.at( i ).second ) ) / 2.0 ).segment( 0, 3 );

//        // Check if observing link end is occulted by body.
//        if( mission_geometry::computeShadowFunction(
//                    linkEndStates.at( linkEndIndices_.at( i ).first ).segment( 0, 3 ), 0.0,
//                    positionOfOccultingBody,
//                    radiusOfOccultingBody_,
//                    linkEndStates.at( linkEndIndices_.at( i ).second ).segment( 0, 3 ) ) < 1.0E-10 )
//        {
            isObservationPossible = true;
            break;
        }
    }//

    return isObservationPossible;
}



}

}
