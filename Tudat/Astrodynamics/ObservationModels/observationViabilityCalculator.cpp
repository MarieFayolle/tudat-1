/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/ObservationModels/observationViabilityCalculator.h"
#include "tudatApplications/master_thesis/Functions/PODmodel.h"

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



//! Function for determining whether the link is within the antenna coverage cone (wrt beamwidth) during the observataion.
bool AntennaCoverageCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                     const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 0;
    Eigen::Vector6d positionOfOccultingBody;

    if ( linkEndIndices_.size() != antennaAngularPositionVector_.size() || linkEndIndices_.size() != antennaBeamwidthVector_.size() ){

        throw std::runtime_error(" Error, unconsistent sizes when comparing link end indices with antenna properties in antenna "
                                 "coverage calculator ");
    }


    else {

        // Iterate over all sets of entries of input vector for which occultation is to be checked.
        for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
        {
            // Get position of occulting body
            positionOfOccultingBody = stateFunctionOfOccultingBody_(
                        ( linkEndTimes.at( linkEndIndices_.at( i ).first ) +
                          linkEndTimes.at( linkEndIndices_.at( i ).second ) ) / 2.0 );


            // Get first antenna inertial orientation
            std::pair< double, double > angularOrientationAntenna = antennaAngularPositionVector_[i].first;

            Eigen::Vector3d sphericalBodyFixedAntennaOrientation = ( Eigen::Vector3d( ) << 1.0,
                                                                     mathematical_constants::PI / 2.0 - angularOrientationAntenna.second,
                                                                     angularOrientationAntenna.first ).finished();

            Eigen::Vector3d cartesianBodyFixedAntennaOrientation = coordinate_conversions::convertSphericalToCartesian(
                        sphericalBodyFixedAntennaOrientation );

            Eigen::Vector3d cartesianInertialAntennaOrientation = rotationalEphemeridesVector_[i].first
                    ->getRotationToBaseFrame( linkEndTimes.at( linkEndIndices_.at( i ).first ) ).toRotationMatrix()
                    * cartesianBodyFixedAntennaOrientation;

            Eigen::Vector3d sphericalInertialAntennaOrientation = coordinate_conversions::convertCartesianToSpherical(
                        cartesianInertialAntennaOrientation );

            std::pair< double, double > firstAntennaInertialOrientation = std::make_pair( sphericalInertialAntennaOrientation[ 2 ],
                    mathematical_constants::PI / 2.0 - sphericalInertialAntennaOrientation[ 1 ]);


            // Get second antenna inertial orientation
            angularOrientationAntenna = antennaAngularPositionVector_[i].second;

            sphericalBodyFixedAntennaOrientation = ( Eigen::Vector3d( ) << 1.0, mathematical_constants::PI / 2.0 - angularOrientationAntenna.second,
                      angularOrientationAntenna.first ).finished();

            cartesianBodyFixedAntennaOrientation = coordinate_conversions::convertSphericalToCartesian( sphericalBodyFixedAntennaOrientation );

            cartesianInertialAntennaOrientation = rotationalEphemeridesVector_[i].first
                    ->getRotationToBaseFrame( linkEndTimes.at( linkEndIndices_.at( i ).first ) ).toRotationMatrix()
                    * cartesianBodyFixedAntennaOrientation;

            sphericalInertialAntennaOrientation = coordinate_conversions::convertCartesianToSpherical( cartesianInertialAntennaOrientation );

            std::pair< double, double > secondAntennaInertialOrientation = std::make_pair( sphericalInertialAntennaOrientation[ 2 ],
                    mathematical_constants::PI / 2.0 - sphericalInertialAntennaOrientation[ 1 ]);



            // Check if observed link end is visible from antenna of the observing link end.
            if( ( visibilityConditionWithBeamwidthConstraintAndSphericalBodyOccultation( antennaBeamwidthVector_[i].first,
                             radiusOfOccultingBody_, firstAntennaInertialOrientation, linkEndStates.at( linkEndIndices_.at( i ).first ),
                             linkEndStates.at( linkEndIndices_.at( i ).second ), positionOfOccultingBody ) == true )

                    && ( visibilityConditionWithBeamwidthConstraintAndSphericalBodyOccultation( antennaBeamwidthVector_[i].second,
                             radiusOfOccultingBody_, secondAntennaInertialOrientation, linkEndStates.at( linkEndIndices_.at( i ).second ),
                             linkEndStates.at( linkEndIndices_.at( i ).first ), positionOfOccultingBody) == true ) ){


                isObservationPossible = 1;
                break;

            }
        }

        return isObservationPossible;

    }

}



}

}
