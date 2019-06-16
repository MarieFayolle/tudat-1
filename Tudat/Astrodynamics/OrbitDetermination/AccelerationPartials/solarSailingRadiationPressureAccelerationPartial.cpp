#include <iostream>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/solarSailingRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for updating partial w.r.t. the bodies' positions.
void SolarSailingRadiationPressurePartial::update( const double currentTime )
{

    if( !( currentTime_ == currentTime ) )
    {
        Eigen::Vector3d unitVectorToSource = radiationPressureAcceleration_->getCurrentVectorToSource( );
        Eigen::Vector3d unitVelocityVector = radiationPressureAcceleration_->getCurrentVelocityVector();
        double distanceToSource = radiationPressureAcceleration_->getCurrentDistanceToSource( );
        double velocityWrtSource = radiationPressureAcceleration_->getCurrentVelocityWrtSource( );
        Eigen::Vector3d currentAcceleration = radiationPressureAcceleration_->getAcceleration( );
        double currentRadiationPressure = radiationPressureAcceleration_->getCurrentRadiationPressure( );
        double currentMass = radiationPressureAcceleration_->getCurrentMass( );
        double currentAccelerationMagnitude = currentAcceleration.norm();

        currentPartialWrtPosition_.setZero( );

        if( currentRadiationPressure > 0.0 && currentAcceleration.norm( ) > 0.0 )
        {
            Eigen::Vector3d unitVectorFromSource = - unitVectorToSource;
            Eigen::Vector3d vectorFromSource = distanceToSource * unitVectorFromSource;
            Eigen::Vector3d velocityVector = velocityWrtSource * unitVelocityVector;

            Eigen::Matrix3d currentSourceUnitVectorPartial =  1.0 / distanceToSource * (
                        Eigen::Matrix3d::Identity( ) - unitVectorFromSource * unitVectorFromSource.transpose( ) );

            Eigen::Matrix3d currentVelocityUnitVectorPartial = 1.0 / velocityWrtSource * (
                        Eigen::Matrix3d::Identity( ) - unitVelocityVector * unitVelocityVector.transpose( ) );          

            double frontEmissivityCoefficient = radiationPressureInterface_->getFrontEmissivityCoefficient();
            double backEmissivityCoefficient = radiationPressureInterface_->getBackEmissivityCoefficient();
            double frontLambertianCoefficient = radiationPressureInterface_->getFrontLambertianCoefficient();
            double backLambertianCoefficient = radiationPressureInterface_->getBackLambertianCoefficient();
            double specularReflectionCoefficient = radiationPressureInterface_->getSpecularReflectionCoefficient();
            double reflectivityCoefficient = radiationPressureInterface_->getReflectivityCoefficient();

            double A = frontLambertianCoefficient * reflectivityCoefficient * ( 1.0 - specularReflectionCoefficient )
                    + ( 1.0 - reflectivityCoefficient ) * ( frontEmissivityCoefficient * frontLambertianCoefficient - backEmissivityCoefficient * backLambertianCoefficient)
                    / ( frontEmissivityCoefficient + backEmissivityCoefficient );

            double cosinusConeAngle = cos( radiationPressureInterface_->getCurrentConeAngle() );
            double sinusConeAngle = sin( radiationPressureInterface_->getCurrentConeAngle() );

            double phi = atan2( ( ( 1.0 - radiationPressureInterface_->getSpecularReflectionCoefficient()
                                    * radiationPressureInterface_->getReflectivityCoefficient() ) * cosinusConeAngle * sinusConeAngle ),
                                ( ( 1.0 + radiationPressureInterface_->getReflectivityCoefficient() * radiationPressureInterface_->getSpecularReflectionCoefficient() )
                                  * std::pow( cosinusConeAngle, 2 ) + A * cosinusConeAngle ) );

            double theta = radiationPressureInterface_->getCurrentConeAngle() - phi;

            // Compute the normalised direction vector of the force defined in the (r,theta,k) frame.
            Eigen::Vector3d forceDirectionLocalFrame = ( Eigen::Vector3d() << cos( theta ),
                                                         sin( theta ) * sin( radiationPressureInterface_->getCurrentClockAngle() ),
                                                         sin( theta ) * cos( radiationPressureInterface_->getCurrentClockAngle() ) ).finished();

//            double forceMagnitude = currentRadiationPressure * radiationPressureAcceleration_->getCurrentArea() * cosinusConeAngle
//                    * std::sqrt( std::pow( cosinusConeAngle * ( 1.0 + radiationPressureInterface_->getReflectivityCoefficient()
//                                                                *  radiationPressureInterface_->getSpecularReflectionCoefficient() ) + A, 2 )
//                                 + std::pow( ( 1.0 - radiationPressureInterface_->getSpecularReflectionCoefficient() * radiationPressureInterface_->getReflectivityCoefficient() ) * sinusConeAngle , 2 ) );

//            double calculatedAccelerationMagnitude = currentRadiationPressure * radiationPressureInterface_->getArea() / currentMass * cosinusConeAngle
//                    * std::sqrt( std::pow( ( 1.0 + specularReflectionCoefficient * reflectivityCoefficient ) * cosinusConeAngle + A , 2 )
//                                 + std::pow( ( 1.0 - specularReflectionCoefficient * reflectivityCoefficient ) * sinusConeAngle, 2 ) );

            Eigen::Matrix< double, 1, 3 > currentRadiationPressurePositionPartial =
                    2.0 * currentRadiationPressure * unitVectorToSource.transpose( ) / ( distanceToSource );


            for ( int j = 0 ; j < 3 ; j++ ){

                Eigen::Matrix3d currentPartialRotationMatrixWrtPosition = Eigen::Matrix3d::Zero();
                Eigen::Matrix3d currentPartialRotationMatrixWrtVelocity = Eigen::Matrix3d::Zero();

                for ( int i = 0 ; i < 3 ; i++ ){

                    currentPartialRotationMatrixWrtPosition( i, 0 ) = currentSourceUnitVectorPartial( j, i );
                    currentPartialRotationMatrixWrtVelocity( i, 0 ) = 0.0;

                    for ( int k = 0 ; k < 3 ; k++ )
                    {
                        if ( k != i )
                        {
                            currentPartialRotationMatrixWrtPosition( i, 1 ) += 2.0 * unitVelocityVector[ i ] *
                                    unitVectorFromSource[ k ] * currentSourceUnitVectorPartial( j, k )
                                    - currentSourceUnitVectorPartial( j, i ) * unitVectorFromSource[ k ] * unitVelocityVector[ k ]
                                    - unitVectorFromSource[ i ] * unitVelocityVector[ k ] * currentSourceUnitVectorPartial( j, k );

                            currentPartialRotationMatrixWrtVelocity( i, 1 ) += currentVelocityUnitVectorPartial( j, i ) * std::pow( unitVectorFromSource[ k ], 2 )
                                    - unitVectorFromSource[ i ] * unitVectorFromSource[ k ] * currentVelocityUnitVectorPartial( j, k );
                        }
                    }

                    currentPartialRotationMatrixWrtPosition( i, 2 ) = currentSourceUnitVectorPartial( j, ( i + 1 ) % 3 ) * unitVelocityVector[ ( i + 2 ) % 3 ]
                            - currentSourceUnitVectorPartial( j, ( i + 2 ) % 3 ) * unitVelocityVector[ ( i + 1 ) % 3 ];

                    currentPartialRotationMatrixWrtVelocity( i, 2 ) = unitVectorFromSource[ ( i + 1 ) % 3 ] * currentVelocityUnitVectorPartial( j, ( i + 2 ) % 3 )
                            - unitVectorFromSource[ ( i + 2 ) % 3 ] * currentVelocityUnitVectorPartial( j, ( i + 1 ) % 3 );

                }


                currentPartialWrtPosition_.block( j, 0, 1, 3 ) = ( currentAccelerationMagnitude
                        * currentPartialRotationMatrixWrtPosition * forceDirectionLocalFrame ).transpose();

                currentPartialWrtVelocity_.block( j, 0, 1, 3 ) = ( currentAccelerationMagnitude
                        * currentPartialRotationMatrixWrtVelocity * forceDirectionLocalFrame ).transpose();

            }

            currentPartialWrtPosition_ += currentAccelerationMagnitude / currentRadiationPressure
                    * currentRadiationPressurePositionPartial.transpose() * ( currentAcceleration.normalized().transpose() );


        }

        currentTime_ = currentTime;

    }
}


}

}

