/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/mutualApproximationPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Compute intermediate variable A, to be used to simplify the expression for the partial w.r.t.
//! link end position of the first and second time derivatives of the right ascension.
double computeIntermediateVariableA( Eigen::Vector3d cartesianPositionVector )
{
    return std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( ) + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
            + cartesianPositionVector.x( );
}

//! Compute partial of intermediate variable A w.r.t. link end position
//! (to be used to simplify the expression for the partial w.r.t. link end position of the first
//! and second time derivatives of the right ascension).
Eigen::Matrix< double, 1, 3 > computePartialOfIntermediateVariableAWrtLinkEndPosition( Eigen::Vector3d cartesianPositionVector )
{
    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    partial( 0 ) = cartesianPositionVector.x( ) / ( std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( )
                                              + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) ) ) + 1.0;
    partial( 1 ) = cartesianPositionVector.y( ) / ( std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( )
                                              + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) ) );
    return partial;
}

//! Compute intermediate variable B, to be used to simplify the expression for the partial w.r.t.
//! link end position of the first and second time derivative of the right ascension.
double computeIntermediateVariableB( Eigen::Vector3d cartesianPositionVector )
{
    double intermediateVariableA = computeIntermediateVariableA( cartesianPositionVector );
    return 1.0 / ( ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
                   * ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) ) );
}

//! Compute partial of intermediate variable B w.r.t. link end position
//! (to be used to simplify the expression for the partial w.r.t. link end position of the first
//! and second time derivatives of the right ascension).
Eigen::Matrix< double, 1, 3 > computePartialOfIntermediateVariableBWrtLinkEndPosition( Eigen::Vector3d cartesianPositionVector )
{
    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    double intermediateVariableA = computeIntermediateVariableA( cartesianPositionVector );
    Eigen::Vector3d partialConstantAWrtPosition = computePartialOfIntermediateVariableAWrtLinkEndPosition( cartesianPositionVector );

    partial( 0 ) = - 4.0 * intermediateVariableA * partialConstantAWrtPosition.x( )
            / ( ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
                * ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
                * ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) ) );

    partial( 1 ) =  - 4.0 * ( intermediateVariableA * partialConstantAWrtPosition.y( ) + cartesianPositionVector.y( ) )
            / ( ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
                * ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) )
                * ( intermediateVariableA * intermediateVariableA + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) ) );

    return partial;
}

//! Compute intermediate variable C, to be used to simplify the expression for the partial w.r.t.
//! link end position of the first and second time derivative of the right ascension and declination.
double computeIntermediateVariableC( Eigen::Vector3d cartesianPositionVector )
{
    return std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( )
                      + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) );
}

//! Compute partial of intermediate variable C w.r.t. link end position
//! (to be used to simplify the expression for the partial w.r.t. link end position of the first
//! and second time derivatives of the right ascension and declination).
Eigen::Matrix< double, 1, 3 > computePartialOfIntermediateVariableCWrtLinkEndPosition( Eigen::Vector3d cartesianPositionVector )
{
    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    partial( 0 ) = cartesianPositionVector.x( ) / std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( )
                                                       + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) );
    partial( 1 ) = cartesianPositionVector.y( ) / std::sqrt( cartesianPositionVector.x( ) * cartesianPositionVector.x( )
                                                       + cartesianPositionVector.y( ) * cartesianPositionVector.y( ) );

    return partial;
}

//! Compute intermediate variable D, to be used to simplify the expression for the partial w.r.t.
//! link end position of the first and second time derivative of the declination.
double computeIntermediateVariableD( Eigen::Vector3d cartesianPositionVector )
{
    return cartesianPositionVector.x( ) * cartesianPositionVector.x( ) + cartesianPositionVector.y( ) * cartesianPositionVector.y( )
            + cartesianPositionVector.z( ) * cartesianPositionVector.z( );
}

//! Compute partial of intermediate variable D w.r.t. link end position
//! (to be used to simplify the expression for the partial w.r.t. link end position of the first
//! and second time derivatives of the declination).
Eigen::Matrix< double, 1, 3 > computePartialOfIntermediateVariableDWrtLinkEndPosition( Eigen::Vector3d cartesianPositionVector )
{
    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    partial( 0 ) = 2.0 * cartesianPositionVector.x( );
    partial( 1 ) = 2.0 * cartesianPositionVector.y( );
    partial( 2 ) = 2.0 * cartesianPositionVector.z( );

    return partial;
}


//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector/*,
        const bool isLinkEndReceiver*/ )
{
//    // Define multiplier of patial vector
//    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - cartesianPositionVector( 1 );
    partial( 1 ) = cartesianPositionVector( 0 );
    partial /= /*( partialMultiplier **/ ( cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) +
                                       cartesianPositionVector( 1 ) * cartesianPositionVector( 1 ) ); // );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d cartesianPositionVector/*,
        const bool isLinkEndReceiver*/ )
{
//    // Define multiplier of patial vector
//    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    double range = cartesianPositionVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - cartesianPositionVector( 0 ) * cartesianPositionVector( 2 );
    partial( 1 ) = - cartesianPositionVector( 1 ) * cartesianPositionVector( 2 );
    partial( 2 ) = cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) + cartesianPositionVector( 1 ) * cartesianPositionVector( 1 );
    partial /= /*partialMultiplier /*/
            ( range * range * std::sqrt( cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) + cartesianPositionVector( 1 ) * cartesianPositionVector( 1 ) ) );

    return partial;
}


double computePartialOfRightAscensionWrtTime( Eigen::Vector3d cartesianPositionVector,
                                              Eigen::Vector3d cartesianVelocityVector )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );
    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );

    double partial = 2.0 * ( std::sqrt( rx * rx + ry * ry ) + rx )
            / ( ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx ) + ry * ry )
            * 1.0 / std::sqrt( rx * rx + ry * ry ) * ( rx * vy - ry * vx );

    return partial;
}


double computePartialOfDeclinationWrtTime( Eigen::Vector3d cartesianPositionVector,
                                           Eigen::Vector3d cartesianVelocityVector )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );
    double rz = cartesianPositionVector.z( );
    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );
    double vz = cartesianVelocityVector.z( );

    double normXandYcomponents = std::sqrt( rx * rx + ry * ry );
    double squaredDistance = rx * rx + ry * ry + rz * rz;

    double partial = ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * vz
            - ( rz / ( std::sqrt( rx * rx + ry * ry ) * squaredDistance ) ) * ( rx * vx + ry * vy + rz * vz );

    return partial;
}


//! Function to compute the partial w.r.t. position of observer or observed object of the first time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector )
{

    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );

    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );

    // Define intermediate variable A = std::sqrt( rx * rx + ry * ry ) + rx and associated partial wrt link end position.
    double intermediateVariableA = computeIntermediateVariableA( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableA = computePartialOfIntermediateVariableAWrtLinkEndPosition( cartesianPositionVector );

    // Define intermediate variable B = 1.0 / ( ( A * A + ry * ry )^2 ).
    double intermediateVariableB = computeIntermediateVariableB( cartesianPositionVector );

    // Define intermediate variable C = std::sqrt( rx * rx + ry * ry ) and associated partial wrt link end position.
    double intermediateVariableC = computeIntermediateVariableC( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableC = computePartialOfIntermediateVariableCWrtLinkEndPosition( cartesianPositionVector );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    double intermediateVariable = intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry );

    partial( 0 ) = 2.0 * ( 1.0 / intermediateVariableC * ( rx * vy - ry * vx ) * intermediateVariableB * ( partialIntermediateVariableA.x( ) * ( intermediateVariableA * intermediateVariableA + ry * ry )
                                                                                   - 2.0 * intermediateVariableA * intermediateVariableA * partialIntermediateVariableA.x( ) )
                           + intermediateVariable * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )
                           + intermediateVariable / intermediateVariableC * vy );


    partial( 1 ) = 2.0 * ( 1.0 / intermediateVariableC * ( rx * vy - ry * vx ) * intermediateVariableB * ( partialIntermediateVariableA.y( ) * ( intermediateVariableA * intermediateVariableA + ry * ry )
                                                                                   - 2.0 * intermediateVariableA * ( intermediateVariableA * partialIntermediateVariableA.y( ) + ry ) )
                           + intermediateVariable * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                           - intermediateVariable / intermediateVariableC * vx );
    return partial;
}


//! Function to compute the partial w.r.t. position of observer or observed object of the first time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );
    double rz = cartesianPositionVector.z( );

    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );
    double vz = cartesianVelocityVector.z( );

    // Define intermediate variable C = std::sqrt( rx * rx + ry * ry ) and associated partial wrt link end position.
    double intermediateVariableC = computeIntermediateVariableC( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableC = computePartialOfIntermediateVariableCWrtLinkEndPosition( cartesianPositionVector );

    // Define intermediate variable D = rx * rx + ry * ry + rz * rz and associated partial wrt link end position.
    double intermediateVariableD = computeIntermediateVariableD( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableE = computePartialOfIntermediateVariableDWrtLinkEndPosition( cartesianPositionVector );


    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    partial( 0 ) = - 1.0 / ( intermediateVariableC * intermediateVariableC ) * vz * partialIntermediateVariableC.x( )
            - ( rx * vx + ry * vy + rz * vz ) * rz * ( - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableE.x( )
                                                       - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( ) )
            - rz / ( intermediateVariableC * intermediateVariableD ) * vx;

    partial( 1 ) = - 1.0 / ( intermediateVariableC * intermediateVariableC ) * vz * partialIntermediateVariableC.y( )
            - ( rx * vx + ry * vy + rz * vz ) * rz * ( - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableE.y( )
                                                       - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( ) )
            - rz / ( intermediateVariableC * intermediateVariableD ) * vy;

    partial( 2 ) = - ( rx * vx + ry * vy + rz * vz ) * ( 1.0 / ( intermediateVariableC * intermediateVariableD )
                                                         - rz / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableE.z( ) )
            - rz / ( intermediateVariableC * intermediateVariableD ) * vz;

    return partial;
}



double computeSecondPartialRightAscensionWrtTime( Eigen::Vector3d cartesianPosition,
                                                  Eigen::Vector3d cartesianVelocity,
                                                  Eigen::Vector3d cartesianAcceleration )
{
    double rx = cartesianPosition[ 0 ];
    double ry = cartesianPosition[ 1 ];
    double vx = cartesianVelocity[ 0 ];
    double vy = cartesianVelocity[ 1 ];
    double ax = cartesianAcceleration[ 0 ];
    double ay = cartesianAcceleration[ 1 ];

    double intermediateVariable = std::sqrt( rx * rx + ry * ry ) + rx;

    double normXandYcomponents = std::sqrt( rx * rx + ry * ry );

//    double partial = 2.0 * ( ( 1.0 / ( ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx ) ) )
//                             * ( 1.0 / sqrt( rx * rx + ry * ry ) ) * ( rx * vy - ry * vx )
//                             * ( ( ry * ry - ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx ) ) * vx
//                                 - 2.0 * ( std::sqrt( rx * rx + ry * ry ) + rx ) * ry * vy
//                                 + ( 1.0 / ( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
//                                 * ( ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx )
//                                     * ( - 2.0 * std::sqrt( rx * rx + ry * ry ) - rx ) - rx * ry * ry ) )
//                             + ( ( std::sqrt( rx * rx + ry * ry ) + rx )
//                                 / ( ( std::sqrt( rx * rx + ry * ry ) + rx ) * ( std::sqrt( rx * rx + ry * ry ) + rx ) + ry * ry ) )
//                             + ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * ay - ry * ax ) );

//    partial = 2.0 * ( 1.0 / ( ( intermediateVariable * intermediateVariable + ry * ry ) * ( intermediateVariable * intermediateVariable + ry * ry ) )
//                      * ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * vy - ry * vx )
//                      * ( ( ry * ry - intermediateVariable * intermediateVariable ) * vx - 2.0 * intermediateVariable * ry * vy
//                      + ( 1.0 / ( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
//                      * ( intermediateVariable * intermediateVariable * ( - 2.0 * std::sqrt( rx * rx + ry * ry ) - rx ) - rx * ry * ry ) )
//                      + intermediateVariable / ( intermediateVariable * intermediateVariable + ry * ry ) * 1.0 / std::sqrt( rx * rx + ry * ry )
//                      * ( rx * ay - ry * ax ) );

    double partial = 2.0 * ( ( 1.0 / ( ( intermediateVariable * intermediateVariable + ry * ry ) * ( intermediateVariable * intermediateVariable + ry * ry ) ) )
                      * ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * vy - ry * vx )
                      * ( ( ry * ry - intermediateVariable * intermediateVariable ) * vx - 2.0 * intermediateVariable * ry * vy
                          + ( 1.0 / ( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
                          * ( intermediateVariable * intermediateVariable * ( - 2.0 * std::sqrt( rx * rx + ry * ry ) - rx ) - rx * ry * ry ) )
                      + intermediateVariable / ( intermediateVariable * intermediateVariable + ry * ry ) * 1.0 / std::sqrt( rx * rx + ry * ry )
                      * ( rx * ay - ry * ax ) );

    return partial;
}


double computeSecondPartialDeclinationWrtTime( Eigen::Vector3d cartesianPosition,
                                              Eigen::Vector3d cartesianVelocity,
                                              Eigen::Vector3d cartesianAcceleration )
{
    double rx = cartesianPosition[ 0 ];
    double ry = cartesianPosition[ 1 ];
    double rz = cartesianPosition[ 2 ];
    double vx = cartesianVelocity[ 0 ];
    double vy = cartesianVelocity[ 1 ];
    double vz = cartesianVelocity[ 2 ];
    double ax = cartesianAcceleration[ 0 ];
    double ay = cartesianAcceleration[ 1 ];
    double az = cartesianAcceleration[ 2 ];

    double squaredDistance = rx * rx + ry * ry + rz * rz;

    double partial = ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * az
            - vz * ( 1.0 / ( rx * rx + ry * ry ) ) * ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
            - ( rx * vx + ry * vy + rz * vz ) * ( 1.0 / ( ( rx * rx + ry * ry ) * squaredDistance * squaredDistance ) )
            * ( std::sqrt( rx * rx + ry * ry ) * squaredDistance * vz
                - rz * ( ( squaredDistance / std::sqrt( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
                         + 2.0 * std::sqrt( rx * rx + ry * ry ) * ( rx * vx + ry * vy + rz * vz ) ) )
            - ( rz / ( std::sqrt( rx * rx + ry * ry ) * squaredDistance ) )
            * ( vx * vx + rx * ax + vy * vy + ry * ay + vz * vz + rz * az );

//    partial = ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * dderivRz
//            - derivRz * ( 1.0 / ( rx * rx + ry * ry ) ) * ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * derivRx + ry * derivRy )
//            - ( rx * derivRx + ry * derivRy + rz * derivRz ) * ( 1.0 / ( rx * rx + ry * ry ) ) * ( 1.0 / ( squaredDistance * squaredDistance ) )
//            * ( std::sqrt( rx * rx + ry * ry ) * squaredDistance * derivRz
//                - rz * ( squaredDistance / std::sqrt( rx * rx + ry * ry ) * ( rx * derivRx + ry * derivRy )
//                         + 2.0 * std::sqrt( rx * rx + ry * ry ) * ( rx * derivRx + ry * derivRy + rz * derivRz ) ) )
//            - rz / std::sqrt( rx * rx + ry * ry ) * 1.0 / squaredDistance
//            * ( derivRx * derivRx + rx * dderivRx + derivRy * derivRy + ry * dderivRy + derivRz * derivRz + rz * dderivRz );

    return partial;
}



//! Function to compute the partial w.r.t. position of observer or observed object of the second time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
        const Eigen::Vector3d& cartesianAccelerationVector,
        const Eigen::Matrix3d& partialCartesianAccelerationWrtPosition,
        const bool wrtTransmitterPosition,
        const bool wrtOtherLinkEnd )
{

    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );

    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );

    double ax = cartesianAccelerationVector.x( );
    double ay = cartesianAccelerationVector.y( );

    Eigen::Vector3d accelerationPartialWrtX = partialCartesianAccelerationWrtPosition.block( 0, 0, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtY = partialCartesianAccelerationWrtPosition.block( 0, 1, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtZ = partialCartesianAccelerationWrtPosition.block( 0, 2, 3, 1 );
    if ( !wrtTransmitterPosition )
    {
        accelerationPartialWrtX *= - 1.0;
        accelerationPartialWrtY *= - 1.0;
        accelerationPartialWrtZ *= - 1.0;
    }

    // Define intermediate variable A = std::sqrt( rx * rx + ry * ry ) + rx and associated partial wrt link end position.
    double intermediateVariableA = computeIntermediateVariableA( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableA = computePartialOfIntermediateVariableAWrtLinkEndPosition( cartesianPositionVector );

    // Define intermediate variable B = 1.0 / ( ( A * A + ry * ry )^2 ) and associated partial wrt link end position.
    double intermediateVariableB = computeIntermediateVariableB( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableB = computePartialOfIntermediateVariableBWrtLinkEndPosition( cartesianPositionVector );

    // Define intermediate variable C = std::sqrt( rx * rx + ry * ry ) and associated partial wrt link end position.
    double intermediateVariableC = computeIntermediateVariableC( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableC = computePartialOfIntermediateVariableCWrtLinkEndPosition( cartesianPositionVector );



    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    double intermediateVariable1 = intermediateVariableB / intermediateVariableC * ( rx * vy - ry * vx );
    double intermediateVariable2 = ( ry * ry - intermediateVariableA * intermediateVariableA ) * vx - 2.0 * intermediateVariableA * ry * vy
            + ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * ( rx * vx + ry * vy )
            * ( intermediateVariableA * intermediateVariableA * ( - 2.0 * intermediateVariableC - rx ) - rx * ry * ry );

    if ( !wrtOtherLinkEnd )
    {
        partial( 0 ) = 2.0 * ( intermediateVariable2
                               * ( partialIntermediateVariableB.x( ) * ( 1.0 / intermediateVariableC ) * ( rx * vy - ry * vx )
                                   + intermediateVariableB * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )
                                   + intermediateVariableB / intermediateVariableC * vy )
                               + intermediateVariable1 * ( - 2.0 * intermediateVariableA * partialIntermediateVariableA.x( ) * vx
                                                           - 2.0 * partialIntermediateVariableA.x( ) * ry * vy
                                                           + ( intermediateVariableA * intermediateVariableA * ( - 2.0 * intermediateVariableC - rx ) - rx * ry * ry )
                                                           * ( ( - 2.0 * partialIntermediateVariableC.x( ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) )
                                                               * ( rx * vx + ry * vy ) + ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * vx )
                                                           + ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * ( rx * vx + ry * vy )
                                                           * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.x( ) * ( - 2.0 * intermediateVariableC - rx )
                                                               + intermediateVariableA * intermediateVariableA * ( - 2.0 * partialIntermediateVariableC.x( ) - 1.0  ) - ry * ry ) )
                               + intermediateVariableB * ( ( intermediateVariableA * intermediateVariableA + ry * ry ) * partialIntermediateVariableA.x( )
                                               - 2.0 * intermediateVariableA * intermediateVariableA * partialIntermediateVariableA.x( ) )
                               * 1.0 / intermediateVariableC * ( rx * ay - ry * ax )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
                               * ( ay + rx * accelerationPartialWrtX.y( ) - ry * accelerationPartialWrtX.x( ) )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * ( rx * ay - ry * ax )
                               * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( ) );


        partial( 1 ) = 2.0 * ( intermediateVariable2 * ( 1.0 / intermediateVariableC * partialIntermediateVariableB.y( ) * ( rx * vy - ry * vx )
                                                         + intermediateVariableB * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                                                         + intermediateVariableB / intermediateVariableC * ( - vx ) )
                               + intermediateVariable1
                               * ( ( 2.0 * ry - 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) ) * vx - 2.0 * intermediateVariableA * vy - 2.0 * ry * partialIntermediateVariableA.y( ) * vy
                                   + ( intermediateVariableA * intermediateVariableA * ( - 2.0 * intermediateVariableC - rx ) - rx * ry * ry )
                                   * ( - 2.0 * partialIntermediateVariableC.y( ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) * ( rx * vx + ry * vy )
                                       + 1.0 / ( intermediateVariableC * intermediateVariableC ) * vy )
                                   + 1.0 / ( intermediateVariableC * intermediateVariableC ) * ( rx * vx + ry * vy )
                                   * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) * ( - 2.0 * intermediateVariableC - rx )
                                       + intermediateVariableA * intermediateVariableA * ( - 2.0 * partialIntermediateVariableC.y( ) ) - 2.0 * rx * ry ) )
                               + intermediateVariableB * ( ( intermediateVariableA * intermediateVariableA + ry * ry ) * partialIntermediateVariableA.y( )
                                               - intermediateVariableA * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) + 2.0 * ry ) )
                               * 1.0 / intermediateVariableC * ( rx * ay - ry * ax )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * ( rx * ay - ry * ax )
                               * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
                               * ( - ax + rx * accelerationPartialWrtY.y( ) - ry * accelerationPartialWrtY.x( ) ) );


        partial( 2 ) = 2.0 * ( intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) / intermediateVariableC
                               * ( rx * accelerationPartialWrtZ.y( ) - ry * accelerationPartialWrtZ.x( ) ) );
    }

    else
    {
        partial( 0 ) = 2.0 * ( /*intermediateVariable2
                               * ( partialIntermediateVariableB.x( ) * ( 1.0 / intermediateVariableC ) * ( rx * vy - ry * vx )
                                   + intermediateVariableB * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )
                                   + intermediateVariableB / intermediateVariableC * vy )
                               + intermediateVariable1 * (*/ /*- 2.0 * intermediateVariableA * partialIntermediateVariableA.x( ) * vx
                                                           - 2.0 * partialIntermediateVariableA.x( ) * ry * vy
                                                           + ( intermediateVariableA * intermediateVariableA * ( - 2.0 * intermediateVariableC - rx ) - rx * ry * ry )
                                                           * ( ( - 2.0 * partialIntermediateVariableC.x( ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) )
                                                               * ( rx * vx + ry * vy ) + ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * vx )
                                                           + ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * ( rx * vx + ry * vy )
                                                           * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.x( ) * ( - 2.0 * intermediateVariableC - rx )
                                                               + intermediateVariableA * intermediateVariableA * ( - 2.0 * partialIntermediateVariableC.x( ) - 1.0  ) - ry * ry ) )*/
//                               + intermediateVariableB * ( ( intermediateVariableA * intermediateVariableA + ry * ry ) * partialIntermediateVariableA.x( )
//                                                           - 2.0 * intermediateVariableA * intermediateVariableA * partialIntermediateVariableA.x( ) )
//                               * 1.0 / intermediateVariableC * ( rx * ay - ry * ax )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
                               * ( ay + rx * accelerationPartialWrtX.y( ) - ry * accelerationPartialWrtX.x( ) )
                               /*+ intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * ( rx * ay - ry * ax )
                               * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )*/ );


        partial( 1 ) = 2.0 * ( /*intermediateVariable2 * ( 1.0 / intermediateVariableC * partialIntermediateVariableB.y( ) * ( rx * vy - ry * vx )
                                                         + intermediateVariableB * ( rx * vy - ry * vx ) * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                                                         + intermediateVariableB / intermediateVariableC * ( - vx ) )
                               + intermediateVariable1
                               * ( ( 2.0 * ry - 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) ) * vx - 2.0 * intermediateVariableA * vy - 2.0 * ry * partialIntermediateVariableA.y( ) * vy
                                   + ( intermediateVariableA * intermediateVariableA * ( - 2.0 * intermediateVariableC - rx ) - rx * ry * ry )
                                   * ( - 2.0 * partialIntermediateVariableC.y( ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) * ( rx * vx + ry * vy )
                                       + 1.0 / ( intermediateVariableC * intermediateVariableC ) * vy )
                                   + 1.0 / ( intermediateVariableC * intermediateVariableC ) * ( rx * vx + ry * vy )
                                   * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) * ( - 2.0 * intermediateVariableC - rx )
                                       + intermediateVariableA * intermediateVariableA * ( - 2.0 * partialIntermediateVariableC.y( ) ) - 2.0 * rx * ry ) )
                               + intermediateVariableB * ( ( intermediateVariableA * intermediateVariableA + ry * ry ) * partialIntermediateVariableA.y( )
                                               - intermediateVariableA * ( 2.0 * intermediateVariableA * partialIntermediateVariableA.y( ) + 2.0 * ry ) )
                               * 1.0 / intermediateVariableC * ( rx * ay - ry * ax )*/
//                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * ( rx * ay - ry * ax )
//                               * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                               + intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
                               * ( - ax + rx * accelerationPartialWrtY.y( ) - ry * accelerationPartialWrtY.x( ) ) );


        partial( 2 ) = 2.0 * ( intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) / intermediateVariableC
                               * ( rx * accelerationPartialWrtZ.y( ) - ry * accelerationPartialWrtZ.x( ) ) );
    }


    if ( !wrtTransmitterPosition )
    {
        partial *= - 1.0;
    }

    return partial;

}



//! Function to compute the partial w.r.t. position of observer or observed object of the second time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
        const Eigen::Vector3d& cartesianAccelerationVector,
        const Eigen::Matrix3d& partialCartesianAccelerationWrtPosition,
        const bool wrtTransmitterPosition,
        const bool wrtOtherLinkEnd )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );
    double rz = cartesianPositionVector.z( );

    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );
    double vz = cartesianVelocityVector.z( );

    double ax = cartesianAccelerationVector.x( );
    double ay = cartesianAccelerationVector.y( );
    double az = cartesianAccelerationVector.z( );

    Eigen::Vector3d accelerationPartialWrtX = partialCartesianAccelerationWrtPosition.block( 0, 0, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtY = partialCartesianAccelerationWrtPosition.block( 0, 1, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtZ = partialCartesianAccelerationWrtPosition.block( 0, 2, 3, 1 );
    if ( !wrtTransmitterPosition )
    {
        accelerationPartialWrtX *= - 1.0;
        accelerationPartialWrtY *= - 1.0;
        accelerationPartialWrtZ *= - 1.0;
    }

    // Define intermediate variable C = std::sqrt( rx * rx + ry * ry ) and associated partial wrt link end position.
    double intermediateVariableC = computeIntermediateVariableC( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableC = computePartialOfIntermediateVariableCWrtLinkEndPosition( cartesianPositionVector );

    // Define intermediate variable D = rx * rx + ry * ry + rz * rz and associated partial wrt link end position.
    double intermediateVariableD = computeIntermediateVariableD( cartesianPositionVector );
    Eigen::Matrix< double, 1, 3 > partialIntermediateVariableD = computePartialOfIntermediateVariableDWrtLinkEndPosition( cartesianPositionVector );


    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    double intermediateVariable1 = intermediateVariableC * intermediateVariableD * vz
            - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * intermediateVariableD + 2.0 * intermediateVariableC * ( rx * vx + ry * vy + rz * vz ) );
    double intermediateVariable2 = vx * vx + rx * ax + vy * vy + ry * ay + vz * vz + rz * az;

    if ( !wrtOtherLinkEnd )
    {
        partial( 0 ) = - 1.0 / ( intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( ) * az + 1.0 / intermediateVariableC * accelerationPartialWrtX.z( )
                - vz / intermediateVariableC * ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * vx
                - vz * ( rx * vx + ry * vy ) * ( - 3.0 ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
                - intermediateVariable1 * ( vx / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                                            + ( rx * vx + ry * vy + rz * vz )
                                            * ( - 2.0 / ( intermediateVariableD * intermediateVariableD * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
                                                - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.x( ) ) )
                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                * ( ( intermediateVariableD * partialIntermediateVariableC.x( ) + intermediateVariableC * partialIntermediateVariableD.x( ) ) * vz
                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.x( )
                             + ( rx * vx + ry * vy ) * intermediateVariableD * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )
                             + intermediateVariableD / intermediateVariableC * vx + 2.0 * partialIntermediateVariableC.x( ) * ( rx * vx + ry * vy + rz * vz )
                             + 2.0 * intermediateVariableC * vx ) )
                - intermediateVariable2 * rz * ( - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
                                                 - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.x( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ax + rx * accelerationPartialWrtX.x( ) + ry * accelerationPartialWrtX.y( )
                                                     + rz * accelerationPartialWrtX.z( ) );


        partial( 1 ) = - 1.0 / ( intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( ) * az + 1.0 / intermediateVariableC * accelerationPartialWrtY.z( )
                - vz / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) * vy
                - vz * ( rx * vx + ry * vy ) * ( - 3.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                - intermediateVariable1 * ( vy / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                                            + ( rx * vx + ry * vy + rz * vz )
                                            * ( - 2.0 / ( intermediateVariableD * intermediateVariableD * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( )
                                                - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.y( ) ) )
                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                * ( ( intermediateVariableD * partialIntermediateVariableC.y( ) + intermediateVariableC * partialIntermediateVariableD.y( ) ) * vz
                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.y( )
                             + ( rx * vx + ry * vy ) * intermediateVariableD * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
                             + intermediateVariableD / intermediateVariableC * vy + 2.0 * partialIntermediateVariableC.y( ) * ( rx * vx + ry * vy + rz * vz )
                             + 2.0 * intermediateVariableC * vy ) )
                - intermediateVariable2 * rz * ( - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( )
                                                 - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.y( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ay + rx * accelerationPartialWrtY.x( ) + ry * accelerationPartialWrtY.y( )
                                                     + rz * accelerationPartialWrtY.z( ) );

        partial( 2 ) = 1.0 / intermediateVariableC * accelerationPartialWrtZ.z( )
                - intermediateVariable1 * ( vz / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                                            + ( rx * vx + ry * vy + rz * vz )
                                            * ( - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.z( ) ) )
                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
                * ( intermediateVariableC * partialIntermediateVariableD.z( ) * vz
                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.z( ) + 2.0 * intermediateVariableC * vz )
                    - ( ( rx * vx + ry * vy ) * intermediateVariableD / intermediateVariableC + 2.0 * intermediateVariableC * ( rx * vx + ry * vy + rz * vz ) ) )
                - intermediateVariable2 / intermediateVariableC * ( 1.0 / intermediateVariableD - rz / ( intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.z( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( az + rx * accelerationPartialWrtZ.x( ) + ry * accelerationPartialWrtZ.y( )
                                                     + rz * accelerationPartialWrtZ.z( ) );

    }

    else
    {
        partial( 0 ) = /*- 1.0 / ( intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( ) * az +*/ 1.0 / intermediateVariableC * accelerationPartialWrtX.z( )
//                - vz / intermediateVariableC * ( 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * vx
//                - vz * ( rx * vx + ry * vy ) * ( - 3.0 ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
//                - intermediateVariable1 * ( vx / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                                            + ( rx * vx + ry * vy + rz * vz )
//                                            * ( - 2.0 / ( intermediateVariableD * intermediateVariableD * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
//                                                - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.x( ) ) )
//                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                * ( ( intermediateVariableD * partialIntermediateVariableC.x( ) + intermediateVariableC * partialIntermediateVariableD.x( ) ) * vz
//                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.x( )
//                             + ( rx * vx + ry * vy ) * intermediateVariableD * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.x( )
//                             + intermediateVariableD / intermediateVariableC * vx + 2.0 * partialIntermediateVariableC.x( ) * ( rx * vx + ry * vy + rz * vz )
//                             + 2.0 * intermediateVariableC * vx ) )
//                - intermediateVariable2 * rz * ( - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.x( )
//                                                 - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.x( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ax + rx * accelerationPartialWrtX.x( ) + ry * accelerationPartialWrtX.y( )
                                                     + rz * accelerationPartialWrtX.z( ) );


        partial( 1 ) = /*- 1.0 / ( intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( ) * az*/ + 1.0 / intermediateVariableC * accelerationPartialWrtY.z( )
//                - vz / ( intermediateVariableC * intermediateVariableC * intermediateVariableC ) * vy
//                - vz * ( rx * vx + ry * vy ) * ( - 3.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
//                - intermediateVariable1 * ( vy / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                                            + ( rx * vx + ry * vy + rz * vz )
//                                            * ( - 2.0 / ( intermediateVariableD * intermediateVariableD * intermediateVariableC * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( )
//                                                - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.y( ) ) )
//                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                * ( ( intermediateVariableD * partialIntermediateVariableC.y( ) + intermediateVariableC * partialIntermediateVariableD.y( ) ) * vz
//                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.y( )
//                             + ( rx * vx + ry * vy ) * intermediateVariableD * ( - 1.0 / ( intermediateVariableC * intermediateVariableC ) ) * partialIntermediateVariableC.y( )
//                             + intermediateVariableD / intermediateVariableC * vy + 2.0 * partialIntermediateVariableC.y( ) * ( rx * vx + ry * vy + rz * vz )
//                             + 2.0 * intermediateVariableC * vy ) )
//                - intermediateVariable2 * rz * ( - 1.0 / ( intermediateVariableD * intermediateVariableC * intermediateVariableC ) * partialIntermediateVariableC.y( )
//                                                 - 1.0 / ( intermediateVariableC * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.y( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ay + rx * accelerationPartialWrtY.x( ) + ry * accelerationPartialWrtY.y( )
                                                     + rz * accelerationPartialWrtY.z( ) );

        partial( 2 ) = 1.0 / intermediateVariableC * accelerationPartialWrtZ.z( )
//                - intermediateVariable1 * ( vz / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                                            + ( rx * vx + ry * vy + rz * vz )
//                                            * ( - 2.0 / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.z( ) ) )
//                - ( rx * vx + ry * vy + rz * vz ) / ( intermediateVariableC * intermediateVariableC * intermediateVariableD * intermediateVariableD )
//                * ( intermediateVariableC * partialIntermediateVariableD.z( ) * vz
//                    - rz * ( ( rx * vx + ry * vy ) / intermediateVariableC * partialIntermediateVariableD.z( ) + 2.0 * intermediateVariableC * vz )
//                    - ( ( rx * vx + ry * vy ) * intermediateVariableD / intermediateVariableC + 2.0 * intermediateVariableC * ( rx * vx + ry * vy + rz * vz ) ) )
//                - intermediateVariable2 / intermediateVariableC * ( 1.0 / intermediateVariableD - rz / ( intermediateVariableD * intermediateVariableD ) * partialIntermediateVariableD.z( ) )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( az + rx * accelerationPartialWrtZ.x( ) + ry * accelerationPartialWrtZ.y( )
                                                     + rz * accelerationPartialWrtZ.z( ) );
    }

    if ( !wrtTransmitterPosition )
    {
        partial *= - 1.0;
    }

    return partial;
}


//! Compute the angle theta used to derive the real solutions of the depressed cubic equation.
double computeAngleThetaRealSolutionsCubicEquation( double intermediateQ,
                                                    double intermediateR )
{
    return std::acos( intermediateR / std::sqrt( - intermediateQ * intermediateQ * intermediateQ ) );
}


////! Compute the partial of the intermediate variable Q (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//computePartialOfIntermediateVariableQWrtLinkEndPosition(
//        Eigen::Vector3d depressedCubicPolynomialCoefficients,
//        std::pair< Eigen::Matrix< double, 3, 3 >, Eigen::Matrix< double, 3, 3 > >
//        partialDepressedPolynomialCoefficientsWrtLinkEndPosition )
//{

//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partials;

//    // Retrieve partials of the second order coefficient of the depressed cubic polynomial, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver).
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//            partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
//            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 0, 0, 1, 3 ),
//                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) );

//    // Retrieve partials of the first order coefficient of the depressed cubic polynomial, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver).
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//            partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
//            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 1, 0, 1, 3 ),
//                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) );

////    // Retrieve partials of the zero order coefficient of the depressed cubic polynomial, for the two link legs
////    // (first transmitter - receiver and second transmitter - receiver).
////    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
////            partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
////            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 2, 0, 1, 3 ),
////                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 2, 0, 1, 3 ) );


//    // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver)
//    partials.first = ( 1.0 / 9.0 ) * ( 3.0 * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first
//                                       - 2.0 * depressedCubicPolynomialCoefficients[ 0 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first );

//    partials.first = ( 1.0 / 9.0 ) * ( 3.0 * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second
//                                       - 2.0 * depressedCubicPolynomialCoefficients[ 0 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second );

//    return partials;
//}


//! Compute the partial of the intermediate variable Q (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
//! The partials are returned by reference.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
void computePartialOfIntermediateVariableQWrtLinkEndPosition(
        Eigen::Vector3d depressedCubicPolynomialCoefficients,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtReceiverPosition )
{
    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 3, 3 > partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtReceiverPosition;
        }


        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );

        // Retrieve partials of the second order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition
                = partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 );

        // Retrieve partials of the first order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition
                = partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition.block( 1, 0, 1, 3 );


        // Compute partials of intermediate variable Q w.r.t. link end position.
        partials = ( 1.0 / 9.0 ) * ( 3.0 * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition
                                           - 2.0 * depressedCubicPolynomialCoefficients[ 0 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfIntermediateVariableQwrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfIntermediateVariableQwrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialOfIntermediateVariableQwrtReceiverPosition = partials;
        }

    }

}



////! Compute the partial of the intermediate variable R (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//computePartialOfIntermediateVariableRWrtLinkEndPosition(
//        Eigen::Vector3d depressedCubicPolynomialCoefficients,
//        std::pair< Eigen::Matrix< double, 3, 3 >, Eigen::Matrix< double, 3, 3 > >
//        partialDepressedPolynomialCoefficientsWrtLinkEndPosition )
//{

//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partials;

//    // Retrieve partials of the second order coefficient of the depressed cubic polynomial, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver).
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//            partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
//            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 0, 0, 1, 3 ),
//                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) );

//    // Retrieve partials of the first order coefficient of the depressed cubic polynomial, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver).
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//            partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
//            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 1, 0, 1, 3 ),
//                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) );

//    // Retrieve partials of the zero order coefficient of the depressed cubic polynomial, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver).
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//            partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
//            std::make_pair( partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 2, 0, 1, 3 ),
//                            partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 2, 0, 1, 3 ) );


//    // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver)
//    partials.first = ( 1.0 / 54.0 ) * ( 9.0 * ( depressedCubicPolynomialCoefficients[ 0 ] * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first
//                                        + depressedCubicPolynomialCoefficients[ 1 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first )
//                                       - 27.0 * partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first
//            - 6.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ]
//            * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.first );

//    partials.second = ( 1.0 / 54.0 ) * ( 9.0 * ( depressedCubicPolynomialCoefficients[ 0 ] * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second
//                                        + depressedCubicPolynomialCoefficients[ 1 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second )
//                                       - 27.0 * partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second
//            - 6.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ]
//            * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition.second );

//    return partials;
//}


//! Compute the partial of the intermediate variable R (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
void computePartialOfIntermediateVariableRWrtLinkEndPosition(
        Eigen::Vector3d depressedCubicPolynomialCoefficients,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtReceiverPosition )
{

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 3, 3 > partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition = partialOfDepressedPolynomialCoefficientsWrtReceiverPosition;
        }


        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );

        // Retrieve partials of the second order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 );

        // Retrieve partials of the first order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition.block( 1, 0, 1, 3 );

        // Retrieve partials of the zero order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndPosition.block( 2, 0, 1, 3 );


        // Compute partials of intermediate variable Q w.r.t. link end position.
        partials = ( 1.0 / 54.0 ) * ( 9.0 * ( depressedCubicPolynomialCoefficients[ 0 ] * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndPosition
                                            + depressedCubicPolynomialCoefficients[ 1 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition )
                                           - 27.0 * partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndPosition
                - 6.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ]
                * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndPosition );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfIntermediateVariableRwrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfIntermediateVariableRwrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialOfIntermediateVariableRwrtReceiverPosition = partials;
        }

    }

}


////! Compute the partial of the angle theta (used to derive the real solutions of the central instant depressed cubic equation) w.r.t. link end position.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//computePartialOfAngleThetaCubicEquationWrtLinkEndPosition(
//        double intermediateVariableQ,
//        double intermediateVariableR,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableQwrtLinkEndPosition,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableRwrtLinkEndPosition )
//{
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partials;


//    // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver)
//    partials.first = std::sqrt( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
//                                / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ + intermediateVariableR * intermediateVariableR ) )
//            * ( 1.0 / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) )
//            * ( std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) * partialsIntermediateVariableRwrtLinkEndPosition.first
//                + 3.0 * intermediateVariableR * intermediateVariableQ * intermediateVariableQ
//                / ( 2.0 * std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) ) * partialsIntermediateVariableQwrtLinkEndPosition.first );

//    partials.second = std::sqrt( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
//                                / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ + intermediateVariableR * intermediateVariableR ) )
//            * ( 1.0 / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) )
//            * ( std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) * partialsIntermediateVariableRwrtLinkEndPosition.second
//                + 3.0 * intermediateVariableR * intermediateVariableQ * intermediateVariableQ
//                / ( 2.0 * std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) ) * partialsIntermediateVariableQwrtLinkEndPosition.second );

//    return partials;
//}


//! Compute the partial of the angle theta (used to derive the real solutions of the central instant depressed cubic equation) w.r.t. link end position.
void computePartialOfAngleThetaCubicEquationWrtLinkEndPosition(
        double intermediateVariableQ,
        double intermediateVariableR,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtReceiverPosition )
{

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndPosition;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtFirstTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtSecondTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtReceiverPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtReceiverPosition;
        }


        // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
        // (first transmitter - receiver and second transmitter - receiver)
        Eigen::Matrix< double, 1, 3 > partials = std::sqrt( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
                                    / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ + intermediateVariableR * intermediateVariableR ) )
                * ( 1.0 / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) )
                * ( std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) * partialsOfIntermediateVariableRwrtLinkEndPosition
                    + 3.0 * intermediateVariableR * intermediateVariableQ * intermediateVariableQ
                    / ( 2.0 * std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) ) * partialsOfIntermediateVariableQwrtLinkEndPosition );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialOfAngleThetaCubicEquationWrtReceiverPosition = partials;
        }

    }

}


////! Compute the partial of the intermediate variable T (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//computePartialOfIntermediateVariableTWrtLinkEndPosition(
//        double intermediateVariableQ,
//        double intermediateVariableR,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableQwrtLinkEndPosition,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableRwrtLinkEndPosition )
//{
//    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
//            + intermediateVariableR * intermediateVariableR;

//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableBetaWrtLinkEndPosition;
//    partialsIntermediateVariableBetaWrtLinkEndPosition.first = 3.0 * intermediateVariableQ * intermediateVariableQ * partialsIntermediateVariableQwrtLinkEndPosition.first
//            + 2.0 * intermediateVariableR * partialsIntermediateVariableRwrtLinkEndPosition.first;
//    partialsIntermediateVariableBetaWrtLinkEndPosition.second = 3.0 * intermediateVariableQ * intermediateVariableQ * partialsIntermediateVariableQwrtLinkEndPosition.second
//            + 2.0 * intermediateVariableR * partialsIntermediateVariableRwrtLinkEndPosition.second;


//    // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver)
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partials;

//    partials.first = 1.0 / ( 3.0 * std::pow( intermediateVariableR - std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
//            * ( partialsIntermediateVariableRwrtLinkEndPosition.first - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
//                * partialsIntermediateVariableBetaWrtLinkEndPosition.first );

//    partials.second = 1.0 / ( 3.0 * std::pow( intermediateVariableR - std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
//            * ( partialsIntermediateVariableRwrtLinkEndPosition.second - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
//                * partialsIntermediateVariableBetaWrtLinkEndPosition.second );

//    return partials;

//}

//! Compute the partial of the intermediate variable T (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
void computePartialOfIntermediateVariableTWrtLinkEndPosition(
        double intermediateVariableQ,
        double intermediateVariableR,
        double intermediateVariableT,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtReceiverPosition )
{
    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndPosition;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtFirstTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtSecondTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtReceiverPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtReceiverPosition;
        }

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableBetaWrtLinkEndPosition =
                3.0 * intermediateVariableQ * intermediateVariableQ * partialsOfIntermediateVariableQwrtLinkEndPosition
                + 2.0 * intermediateVariableR * partialsOfIntermediateVariableRwrtLinkEndPosition;


        // Compute partials of intermediate variable Q w.r.t. link end position.
        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
        if ( intermediateVariableT >= 0.0 )
        {
            partials = 1.0 / ( 3.0 * std::pow( intermediateVariableR - std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndPosition - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndPosition );
        }
        else
        {
            partials = 1.0 / ( 3.0 * std::pow( - ( intermediateVariableR - std::sqrt( intermediateVariableBeta ) ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndPosition - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndPosition );
        }


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfIntermediateVariableTwrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfIntermediateVariableTwrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfIntermediateVariableTwrtReceiverPosition = partials;
        }
    }

}



////! Compute the partial of the intermediate variable S (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
//std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > >
//computePartialOfIntermediateVariableSWrtLinkEndPosition(
//        double intermediateVariableQ,
//        double intermediateVariableR,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableQwrtLinkEndPosition,
//        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableRwrtLinkEndPosition )
//{
//    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
//            + intermediateVariableR * intermediateVariableR;

//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableBetaWrtLinkEndPosition;
//    partialsIntermediateVariableBetaWrtLinkEndPosition.first = 3.0 * intermediateVariableQ * intermediateVariableQ * partialsIntermediateVariableQwrtLinkEndPosition.first
//            + 2.0 * intermediateVariableR * partialsIntermediateVariableRwrtLinkEndPosition.first;
//    partialsIntermediateVariableBetaWrtLinkEndPosition.second = 3.0 * intermediateVariableQ * intermediateVariableQ * partialsIntermediateVariableQwrtLinkEndPosition.second
//            + 2.0 * intermediateVariableR * partialsIntermediateVariableRwrtLinkEndPosition.second;


//    // Compute partials of intermediate variable Q w.r.t. link end position, for the two link legs
//    // (first transmitter - receiver and second transmitter - receiver)
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partials;

//    partials.first = 1.0 / ( 3.0 * std::pow( intermediateVariableR + std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
//            * ( partialsIntermediateVariableRwrtLinkEndPosition.first + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
//                * partialsIntermediateVariableBetaWrtLinkEndPosition.first );

//    partials.second = 1.0 / ( 3.0 * std::pow( intermediateVariableR + std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
//            * ( partialsIntermediateVariableRwrtLinkEndPosition.second + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
//                * partialsIntermediateVariableBetaWrtLinkEndPosition.second );

//    return partials;

//}


//! Compute the partial of the intermediate variable S (used to find the roots of the central instant depressed cubic polynomial) w.r.t. link end position.
void computePartialOfIntermediateVariableSWrtLinkEndPosition(
        double intermediateVariableQ,
        double intermediateVariableR,
        double intermediateVariableS,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtFirstTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtSecondTransmitterPosition,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtReceiverPosition )
{
    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndPosition;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtFirstTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtSecondTransmitterPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfIntermediateVariableQwrtLinkEndPosition = partialsOfIntermediateVariableQwrtReceiverPosition;
            partialsOfIntermediateVariableRwrtLinkEndPosition = partialsOfIntermediateVariableRwrtReceiverPosition;
        }

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableBetaWrtLinkEndPosition = 3.0 * intermediateVariableQ * intermediateVariableQ
                * partialsOfIntermediateVariableQwrtLinkEndPosition + 2.0 * intermediateVariableR * partialsOfIntermediateVariableRwrtLinkEndPosition;

        // Compute partials of intermediate variable Q w.r.t. link end position.
        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
        if ( intermediateVariableS >= 0.0 )
        {
            partials = 1.0 / ( 3.0 * std::pow( intermediateVariableR + std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndPosition + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndPosition );
        }
        else
        {
            partials = 1.0 / ( 3.0 * std::pow( - ( intermediateVariableR + std::sqrt( intermediateVariableBeta ) ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndPosition + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndPosition );
        }


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfIntermediateVariableSwrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfIntermediateVariableSwrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfIntermediateVariableSwrtReceiverPosition = partials;
        }
    }

}


Eigen::Vector2d MutualApproximationScaling::computeRelativePositionInInstrumentalFrame( )
{
    Eigen::Vector2d relativePositionInInstrumentalFrame = Eigen::Vector2d::Zero( );

    relativePositionInInstrumentalFrame[ 0 ] = ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
            * std::cos( averageDeclination_ );
    relativePositionInInstrumentalFrame[ 1 ] = declinationSecondTransmitter_ - declinationFirstTransmitter_;

    return relativePositionInInstrumentalFrame;
}


Eigen::Vector2d MutualApproximationScaling::computeRelativeVelocityInInstrumentalFrame( )
{
//    double averageDeclination = ( declinationFirstObject + declinationSecondObject ) / 2.0;

//    double partialOfRightAscensionFirstTransmitter = computePartialOfRightAscensionWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

//    double partialOfRightAscensionSecondTransmitter = computePartialOfRightAscensionWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ) );

//    double partialOfDeclinationFirstTransmitter = computePartialOfDeclinationWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

//    double partialOfDeclinationSecondTransmitter = computePartialOfDeclinationWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ) );

    Eigen::Vector2d relativeVelocityInInstrumentalFrame = Eigen::Vector2d::Zero( );

    relativeVelocityInInstrumentalFrame[ 0 ] = ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( averageDeclination_ ) * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    relativeVelocityInInstrumentalFrame[ 1 ] = partialOfDeclinationSecondTransmitterWrtTime_ - partialOfDeclinationFirstTransmitterWrtTime_;

    return relativeVelocityInInstrumentalFrame;
}


Eigen::Vector2d MutualApproximationScaling::computeRelativeAccelerationInInstrumentalFrame(
        Eigen::Vector3d cartesianAccelerationFirstTransmitterWrtReceiver,
        Eigen::Vector3d cartesianAccelerationSecondTransmitterWrtReceiver )
{

//    double partialOfRightAscensionFirstTransmitterWrtTime = computePartialOfRightAscensionWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

//    double partialOfRightAscensionSecondTransmitterWrtTime = computePartialOfRightAscensionWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ) );

//    double partialOfDeclinationFirstTransmitterWrtTime = computePartialOfDeclinationWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

//    double partialOfDeclinationSecondTransmitterWrtTime = computePartialOfDeclinationWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ) );

//    double secondPartialOfRightAscensionFirstTransmitterWrtTime = computeSecondPartialRightAscensionWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationFirstTransmitterWrtReceiver );

//    double secondPartialOfRightAscensionSecondTransmitterWrtTime = computeSecondPartialRightAscensionWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                cartesianAccelerationSecondTransmitterWrtReceiver );

//    double secondPartialOfDeclinationFirstTransmitterWrtTime = computeSecondPartialDeclinationWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationFirstTransmitterWrtReceiver );

//    double secondPartialOfDeclinationSecondTransmitterWrtTime = computeSecondPartialDeclinationWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                cartesianAccelerationSecondTransmitterWrtReceiver );

    Eigen::Vector2d relativeAccelerationInInstrumentalFrame = Eigen::Vector2d::Zero( );

    relativeAccelerationInInstrumentalFrame[ 0 ] =
            ( secondPartialOfRightAscensionSecondTransmitterWrtTime_ - secondPartialOfRightAscensionFirstTransmitterWrtTime_ ) * std::cos( averageDeclination_ )
            - ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ ) * std::sin( averageDeclination_ )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            - ( ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ ) * std::cos( averageDeclination_ )
            - ( ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 ) * std::sin( averageDeclination_ )
            * ( secondPartialOfDeclinationFirstTransmitterWrtTime_ + secondPartialOfDeclinationSecondTransmitterWrtTime_ );

    relativeAccelerationInInstrumentalFrame[ 1 ] =
            secondPartialOfDeclinationSecondTransmitterWrtTime_ - secondPartialOfDeclinationFirstTransmitterWrtTime_;

    return relativeAccelerationInInstrumentalFrame;
}



//std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > MutualApproximationScaling::computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPosition( )
//{

//    std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > partialsRelativePositionWrtLinkEndPosition
//            = std::make_pair( Eigen::Matrix< double, 2, 3 >::Zero( ), Eigen::Matrix< double, 2, 3 >::Zero( ) );


//    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
//    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );


//    /// First transmitter - receiver link.

//    // Compute partials of right ascension and declination of the first transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

//    // Compute partials of the relative position (X and Y coordinates in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - first transmitter link.
//    partialsRelativePositionWrtLinkEndPosition.first.block( 0, 0, 1, 3 ) =
//            - partialOfRightAscensionWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfDeclinationWrtPositionFirstTransmitter;

//    partialsRelativePositionWrtLinkEndPosition.first.block( 1, 0, 1, 3 ) = - partialOfDeclinationWrtPositionFirstTransmitter;



//    /// Second transmitter - receiver link.

//    // Compute partials of right ascension and declination of the second transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

//    // Compute partials of the relative position (X and Y in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - second transmitter link.
//    partialsRelativePositionWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) =
//            partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfDeclinationWrtPositionSecondTransmitter;

//    partialsRelativePositionWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) = partialOfDeclinationWrtPositionSecondTransmitter;


//    return partialsRelativePositionWrtLinkEndPosition;
//}



void MutualApproximationScaling::computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPosition( )
{

    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );


    /// First transmitter - receiver link.

    // Compute partials of right ascension and declination of the first transmitter w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

    // Compute partials of the relative position (X and Y coordinates in the instrumental frame of the receiver) w.r.t. first transmitter position.
    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_.block
    /*partialsRelativePositionWrtLinkEndPosition.first.block*/( 0, 0, 1, 3 ) =
            - partialOfRightAscensionWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfDeclinationWrtPositionFirstTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_
    /*partialsRelativePositionWrtLinkEndPosition.first*/.block( 1, 0, 1, 3 ) = - partialOfDeclinationWrtPositionFirstTransmitter;

    // Compute contribution of the link first transmitter-receiver to the partials of relative position
    // (X and Y coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_ = - partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;



    /// Second transmitter - receiver link.

    // Compute partials of right ascension and declination of the second transmitter w.r.t. link end position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    // Compute partials of the relative position (X and Y in the instrumental frame of the receiver) w.r.t. second transmitter position.
    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_
    /*partialsRelativePositionWrtLinkEndPosition.second*/.block( 0, 0, 1, 3 ) =
            partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_
    /*partialsRelativePositionWrtLinkEndPosition.second*/.block( 1, 0, 1, 3 ) = partialOfDeclinationWrtPositionSecondTransmitter;

    // Compute contribution of the link second transmitter-receiver to the partials of relative position
    // (X and Y coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_ -= partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;


}




//std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > MutualApproximationScaling::computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPosition( )
//{
//    std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > partialsRelativeVelocityWrtLinkEndPosition
//            = std::make_pair( Eigen::Matrix< double, 2, 3 >::Zero( ), Eigen::Matrix< double, 2, 3 >::Zero( ) );


//    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
//    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

//    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
//    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


//    // Compute partials of right ascension and declination of the first transmitter w.r.t. time.
//    double partialRightAscensionWrtTimeFirstTransmitter = computePartialOfRightAscensionWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );
//    double partialDeclinationWrtTimeFirstTransmitter = computePartialOfDeclinationWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

//    // Compute partials of right ascension and declination of the second transmitter w.r.t. time.
//    double partialRightAscensionWrtTimeSecondTransmitter = computePartialOfRightAscensionWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );
//    double partialDeclinationWrtTimeSecondTransmitter = computePartialOfDeclinationWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


//    /// First transmitter - receiver link.

//    // Compute partials of right ascension and declination of the first transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );


//    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

//    // Compute partial of first time derivative of the declination w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );


//    // Compute partials of the relative velocity (Vx and Vy coordinates in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - first transmitter link.

//    partialsRelativeVelocityWrtLinkEndPosition.first.block( 0, 0, 1, 3 ) =
//            - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
//            - 1.0 / 2.0 * ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter )
//            * partialOfDeclinationWrtPositionFirstTransmitter * std::sin( averageDeclination_ )
//            + 1.0 / 2.0 * partialOfRightAscensionWrtPositionFirstTransmitter * std::sin( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
//            * partialOfDeclinationWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

//    partialsRelativeVelocityWrtLinkEndPosition.first.block( 1, 0, 1, 3 ) = - partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;



//    /// Second transmitter - receiver link.

//    // Compute partials of right ascension and declination of the second transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );


//    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the second transmitter.
//    partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );

//    // Compute partial of first time derivative of the declination w.r.t. link end position, for the second transmitter.
//    partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


//    // Compute partials of the relative velocity (Vx and Vy coordinates in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - second transmitter link.

//    partialsRelativeVelocityWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) =
//            partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
//            - 1.0 / 2.0 * ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter )
//            * partialOfDeclinationWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
//            - 1.0 / 2.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
//            * partialOfDeclinationWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

//    partialsRelativeVelocityWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) = partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;


//    return partialsRelativeVelocityWrtLinkEndPosition;
//}


void MutualApproximationScaling::computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPosition( )
{

    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


    // Compute partials of right ascension and declination of the first transmitter w.r.t. time.
    double partialRightAscensionWrtTimeFirstTransmitter = computePartialOfRightAscensionWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );
    double partialDeclinationWrtTimeFirstTransmitter = computePartialOfDeclinationWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

    // Compute partials of right ascension and declination of the second transmitter w.r.t. time.
    double partialRightAscensionWrtTimeSecondTransmitter = computePartialOfRightAscensionWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );
    double partialDeclinationWrtTimeSecondTransmitter = computePartialOfDeclinationWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


    /// First transmitter - receiver link.

    // Compute partials of right ascension and declination of the first transmitter w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );


    // Compute partial of first time derivative of the right ascension w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

    // Compute partial of first time derivative of the declination w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );


    // Compute partials of the relative velocity (Vx and Vy coordinates in the instrumental frame of the receiver) w.r.t. first transmitter position.

    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_
    /*partialsRelativeVelocityWrtLinkEndPosition.first*/.block( 0, 0, 1, 3 ) =
            - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter )
            * partialOfDeclinationWrtPositionFirstTransmitter * std::sin( averageDeclination_ )
            + 1.0 / 2.0 * partialOfRightAscensionWrtPositionFirstTransmitter * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
            * partialOfDeclinationWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_
    /*partialsRelativeVelocityWrtLinkEndPosition.first*/.block( 1, 0, 1, 3 ) = - partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

    // Compute contribution of the link first transmitter-receiver to the partials of relative velocity
    // (Vx and Vy coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_ = - partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_;



    /// Second transmitter - receiver link.

    // Compute partials of right ascension and declination of the second transmitter w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );


    // Compute partial of first time derivative of the right ascension w.r.t. second transmitter position.
    partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );

    // Compute partial of first time derivative of the declination w.r.t. second transmitter position.
    partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


    // Compute partials of the relative velocity (Vx and Vy coordinates in the instrumental frame of the receiver) w.r.t. second transmitter position.

    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_
    /*partialsRelativeVelocityWrtLinkEndPosition.second*/.block( 0, 0, 1, 3 ) =
            partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter )
            * partialOfDeclinationWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
            - 1.0 / 2.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
            * partialOfDeclinationWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_
    /*partialsRelativeVelocityWrtLinkEndPosition.second*/.block( 1, 0, 1, 3 ) = partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

    // Compute contribution of the link second transmitter-receiver to the partials of relative velocity
    // (Vx and Vy coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_ -= partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;

}



//std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > MutualApproximationScaling::computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
//        Eigen::Vector3d relativeAccelerationFirstTransmitterWrtReceiver,
//        Eigen::Vector3d relativeAccelerationSecondTransmitterWrtReceiver,
//        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtLinkEndPosition,
//        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtLinkEndPosition )
//{
//    std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > partialsRelativeAccelerationWrtLinkEndPosition
//            = std::make_pair( Eigen::Matrix< double, 2, 3 >::Zero( ), Eigen::Matrix< double, 2, 3 >::Zero( ) );


//    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
//    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

//    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
//    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


//    // Compute partials of right ascension and declination of the first transmitter w.r.t. time.
//    double partialRightAscensionWrtTimeFirstTransmitter = computePartialOfRightAscensionWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );
//    double partialDeclinationWrtTimeFirstTransmitter = computePartialOfDeclinationWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

//    // Compute partials of right ascension and declination of the second transmitter w.r.t. time.
//    double partialRightAscensionWrtTimeSecondTransmitter = computePartialOfRightAscensionWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );
//    double partialDeclinationWrtTimeSecondTransmitter = computePartialOfDeclinationWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


//    // Compute second partials of right ascension and declination of the first transmitter w.r.t. time.
//    double secondPartialRightAscensionWrtTimeFirstTransmitter = computeSecondPartialRightAscensionWrtTime(
//                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver, relativeAccelerationFirstTransmitterWrtReceiver );
//    double secondPartialDeclinationWrtTimeFirstTransmitter = computeSecondPartialDeclinationWrtTime(
//                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver, relativeAccelerationFirstTransmitterWrtReceiver );

//    // Compute second partials of right ascension and declination of the second transmitter w.r.t. time.
//    double secondPartialRightAscensionWrtTimeSecondTransmitter = computeSecondPartialRightAscensionWrtTime(
//                relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver, relativeAccelerationSecondTransmitterWrtReceiver );
//    double secondPartialDeclinationWrtTimeSecondTransmitter = computeSecondPartialDeclinationWrtTime(
//                relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver, relativeAccelerationSecondTransmitterWrtReceiver );




//    /// First transmitter - receiver link.

//    // Compute partials of right ascension and declination of the first transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );



//    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

//    // Compute partial of first time derivative of the declination w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );



//    // Compute partial of second time derivative of the right ascension w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
//                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtLinkEndPosition );

//    // Compute partial of second time derivative of the declination w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
//                                                                               relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtLinkEndPosition );



//    // Compute partials of the relative acceleration (Ax and Ay coordinates in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - first transmitter link.

//    partialsRelativeAccelerationWrtLinkEndPosition.first.block( 0, 0, 1, 3 ) =
//            - partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
//            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
//            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionFirstTransmitter

//            - ( - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition ) * std::sin( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * 1.0 / 2.0 * partialOfDeclinationWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
//            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

//            - 1.0 / 4.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::cos( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
//            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * partialOfDeclinationWrtPositionFirstTransmitter * std::sin( averageDeclination_ )

//            - 1.0 / 2.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::sin( averageDeclination_ )
//            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionFirstTransmitter
//            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition ;

//    partialsRelativeAccelerationWrtLinkEndPosition.first.block( 1, 0, 1, 3 ) = - partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;




//    /// Second transmitter - receiver link.

//    // Compute partials of right ascension and declination of the second transmitter w.r.t. link end position.
//    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
//            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

//    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
//            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );



//    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the second transmitter.
//    partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );

//    // Compute partial of first time derivative of the declination w.r.t. link end position, for the second transmitter.
//    partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );



//    // Compute partial of second time derivative of the right ascension w.r.t. link end position, for the second transmitter.
//    partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
//                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtLinkEndPosition );

//    // Compute partial of second time derivative of the declination w.r.t. link end position, for the second transmitter.
//    partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
//                                                                               relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtLinkEndPosition );



//    // Compute partials of the relative acceleration (Ax and Ay coordinates in the instrumental frame of the receiver)
//    // w.r.t. link end position, for the receiver - second transmitter link.

//    partialsRelativeAccelerationWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) =
//            partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition * std::cos( averageDeclination_ )
//            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
//            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionSecondTransmitter

//            - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::sin( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * 1.0 / 2.0 * partialOfDeclinationWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
//            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

//            - 1.0 / 4.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
//            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
//            * partialOfDeclinationWrtPositionSecondTransmitter * std::sin( averageDeclination_ )

//            - 1.0 / 2.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
//            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionSecondTransmitter
//            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
//            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
//            * partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;


//    partialsRelativeAccelerationWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) = partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;



//    return partialsRelativeAccelerationWrtLinkEndPosition;
//}



void MutualApproximationScaling::computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
        Eigen::Vector3d relativeAccelerationFirstTransmitterWrtReceiver,
        Eigen::Vector3d relativeAccelerationSecondTransmitterWrtReceiver,
//            Eigen::Matrix3d partialAccelerationFirstTransmitterWrtLinkEndPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterPosition,
//            Eigen::Matrix3d partialAccelerationSecondTransmitterWrtLinkEndPosition,
        Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition )
{
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition = Eigen::MatrixXd::Zero( 2, 3 );
    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition = Eigen::MatrixXd::Zero( 2, 3 );
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition = Eigen::MatrixXd::Zero( 2, 3 );

//    std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > partialsRelativeAccelerationWrtLinkEndPosition
//            = std::make_pair( Eigen::Matrix< double, 2, 3 >::Zero( ), Eigen::Matrix< double, 2, 3 >::Zero( ) );


    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


    // Compute partials of right ascension and declination of the first transmitter w.r.t. time.
    double partialRightAscensionWrtTimeFirstTransmitter =
            computePartialOfRightAscensionWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );
    double partialDeclinationWrtTimeFirstTransmitter =
            computePartialOfDeclinationWrtTime( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

    // Compute partials of right ascension and declination of the second transmitter w.r.t. time.
    double partialRightAscensionWrtTimeSecondTransmitter =
            computePartialOfRightAscensionWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );
    double partialDeclinationWrtTimeSecondTransmitter =
            computePartialOfDeclinationWrtTime( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );


    // Compute second partials of right ascension and declination of the first transmitter w.r.t. time.
    double secondPartialRightAscensionWrtTimeFirstTransmitter = computeSecondPartialRightAscensionWrtTime(
                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver, relativeAccelerationFirstTransmitterWrtReceiver );
    double secondPartialDeclinationWrtTimeFirstTransmitter = computeSecondPartialDeclinationWrtTime(
                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver, relativeAccelerationFirstTransmitterWrtReceiver );

    // Compute second partials of right ascension and declination of the second transmitter w.r.t. time.
    double secondPartialRightAscensionWrtTimeSecondTransmitter = computeSecondPartialRightAscensionWrtTime(
                relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver, relativeAccelerationSecondTransmitterWrtReceiver );
    double secondPartialDeclinationWrtTimeSecondTransmitter = computeSecondPartialDeclinationWrtTime(
                relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver, relativeAccelerationSecondTransmitterWrtReceiver );




    /// First transmitter - receiver link.

    // Compute partials of right ascension and declination of the first transmitter w.r.t. link end position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );



    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the first transmitter.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );

    // Compute partial of first time derivative of the declination w.r.t. link end position, for the first transmitter.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver );



//    // Compute partial of second time derivative of the right ascension w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
//                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtLinkEndPosition );

    /// do it separately for partials w.r.t. receiver and transmitter.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtTransmitterPosition, true );
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtReceiverPosition, false );

    // Multiply by factor -1 because this partial specifically is w.r.t. receiver state while the others are all w.r.t. transmitter state.
    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter state
    partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition *= - 1.0;


    // Compute partials  of the right ascension of the second transmitter w.r.t. first transmitter position
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition, true, true );

//    // Compute partial of second time derivative of the declination w.r.t. link end position, for the first transmitter.
//    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
//                                                                               relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtLinkEndPosition );

    /// do it separately for partials w.r.t. receiver and transmitter.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtTransmitterPosition, true );
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtReceiverPosition, false );

    // Multiply by factor -1 because this partial specifically is w.r.t. receiver state while the others are all w.r.t. transmitter state.
    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter state
    partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition *= - 1.0;


    // Compute partials  of the declination of the second transmitter w.r.t. first transmitter position
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                               relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition, true, true );



    // Compute partials of the relative acceleration (Ax and Ay coordinates in the instrumental frame of the receiver)
    // w.r.t. link end position, for the receiver - first transmitter link.

//    partialsRelativeAccelerationWrtLinkEndPosition.first.block( 0, 0, 1, 3 ) =
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition.block( 0, 0, 1, 3 ) =
            ( partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition - /*partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition*/ partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition ) * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionFirstTransmitter

            - ( - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition ) * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * 1.0 / 2.0 * partialOfDeclinationWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

            - 1.0 / 4.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfDeclinationWrtPositionFirstTransmitter * std::sin( averageDeclination_ )

            - 1.0 / 2.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::sin( averageDeclination_ )
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionFirstTransmitter
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * /*partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition*/ ( partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition
                                                                               + partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition) ;

//    partialsRelativeAccelerationWrtLinkEndPosition.first.block( 1, 0, 1, 3 ) = - partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition.block( 1, 0, 1, 3 ) = partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition
            - partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition /*partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition*/;



    /// SAME W.R.T. RECEIVER POSITION FOR FIRST TRANSMITTER - RECEIVER LINK
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition.block( 0, 0, 1, 3 ) = - 1.0 * (
            - partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionFirstTransmitter

            - ( - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition ) * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * 1.0 / 2.0 * partialOfDeclinationWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

            - 1.0 / 4.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfDeclinationWrtPositionFirstTransmitter * std::sin( averageDeclination_ )

            - 1.0 / 2.0 * ( - partialOfRightAscensionWrtPositionFirstTransmitter ) * std::sin( averageDeclination_ )
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionFirstTransmitter
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition ) ;

//    partialsRelativeAccelerationWrtLinkEndPosition.first.block( 1, 0, 1, 3 ) = - partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition.block( 1, 0, 1, 3 ) = - 1.0 * ( - partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition );




    /// Second transmitter - receiver link.

    // Compute partials of right ascension and declination of the second transmitter w.r.t. link end position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );



    // Compute partial of first time derivative of the right ascension w.r.t. link end position, for the second transmitter.
    partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );

    // Compute partial of first time derivative of the declination w.r.t. link end position, for the second transmitter.
    partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver );



//    // Compute partial of second time derivative of the right ascension w.r.t. link end position, for the second transmitter.
//    partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
//                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtLinkEndPosition );

    /// do it separately for partials w.r.t. receiver and transmitter.
    partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtTransmitterPosition, true );
    partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtReceiverPosition, false );

    // Multiply by factor -1 because this partial specifically is w.r.t. receiver state while the others are all w.r.t. transmitter state.
    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter state
    partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition *= - 1.0;

    // Compute partials  of the right ascension of the first transmitter w.r.t. second transmitter position
    partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );


//    // Compute partial of second time derivative of the declination w.r.t. link end position, for the second transmitter.
//    partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition =
//            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
//                                                                               relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtLinkEndPosition );

    /// do it separately for partials w.r.t. receiver and transmitter.
    partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtTransmitterPosition, true );
    partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver, partialAccelerationSecondTransmitterWrtReceiverPosition, false );

    // Multiply by factor -1 because this partial specifically is w.r.t. receiver state while the others are all w.r.t. transmitter state.
    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter state
    partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition *= - 1.0;

    // Compute partials  of the declination of the first transmitter w.r.t. second transmitter position
    partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                               relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );



    // Compute partials of the relative acceleration (Ax and Ay coordinates in the instrumental frame of the receiver)
    // w.r.t. link end position, for the receiver - second transmitter link.

//    partialsRelativeAccelerationWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) =
    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition.block( 0, 0, 1, 3 ) =
            ( /*partialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition*/ partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition - partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition ) * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionSecondTransmitter

            - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * 1.0 / 2.0 * partialOfDeclinationWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

            - 1.0 / 4.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfDeclinationWrtPositionSecondTransmitter * std::sin( averageDeclination_ )

            - 1.0 / 2.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionSecondTransmitter
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * /*partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition*/ ( partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition + partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition );

//    partialsRelativeAccelerationWrtLinkEndPosition.second.block( 1, 0, 1, 3 ) = partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition;
    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition.block( 1, 0, 1, 3 ) = partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition - partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition /*partialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition*/;


    /// SAME W.R.T. RECEIVER POSITION FOR SECOND TRANSMITTER - RECEIVER LINK (ADDED TO CONTRIBUTION OF THE FIRST TRANSMITTER - RECEIVER LINK)
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition.block( 0, 0, 1, 3 ) += - 1.0 * (
            partialOfSecondTimeDerivativeRightAscensionWrtReceiverPosition * std::cos( averageDeclination_ )
            - 1.0 / 2.0 * ( secondPartialRightAscensionWrtTimeSecondTransmitter - secondPartialRightAscensionWrtTimeFirstTransmitter )
            * std::sin( averageDeclination_ ) * partialOfDeclinationWrtPositionSecondTransmitter

            - partialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition * std::sin( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * 1.0 / 2.0 * partialOfDeclinationWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            - ( partialRightAscensionWrtTimeSecondTransmitter - partialRightAscensionWrtTimeFirstTransmitter ) * std::sin( averageDeclination_ )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition

            - 1.0 / 4.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter ) * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition * std::cos( averageDeclination_ )
            + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 8.0  * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * ( partialDeclinationWrtTimeFirstTransmitter + partialDeclinationWrtTimeSecondTransmitter )
            * partialOfDeclinationWrtPositionSecondTransmitter * std::sin( averageDeclination_ )

            - 1.0 / 2.0 * partialOfRightAscensionWrtPositionSecondTransmitter * std::sin( averageDeclination_ )
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0 * partialOfDeclinationWrtPositionSecondTransmitter
            * ( secondPartialDeclinationWrtTimeFirstTransmitter + secondPartialDeclinationWrtTimeSecondTransmitter ) * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition );

    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition.block( 1, 0, 1, 3 ) += - 1.0 * ( partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition );



//    return partialsRelativeAccelerationWrtLinkEndPosition;
}



//std::pair< Eigen::Matrix< double, 4, 3 >, Eigen::Matrix< double, 4, 3 > >
//MutualApproximationScaling::computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition( )
//{
//    // Set partial matrix.
//    std::pair< Eigen::Matrix< double, 4, 3 >, Eigen::Matrix< double, 4, 3 > > partials
//            = std::make_pair( Eigen::Matrix< double, 4, 3 >::Zero( ), Eigen::Matrix< double, 4, 3 >::Zero( ) );

//    // First consider the first transmitter - receiver leg independently, then the second transmitter - receiver leg.
//    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition;
//    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition;
//    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition;


//    for ( unsigned int currentLegIndex = 0 ; currentLegIndex < 2 ; currentLegIndex++ )
//    {

//        if ( currentLegIndex == 0 ) // If first transmitter - receiver leg
//        {
//            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition_.first;
//            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_.first;
//            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_.first;
//        }
//        else // If second transmitter - receiver leg
//        {
//            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition_.second;
//            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_.second;
//            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition =
//                    partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_.second;
//        }


//        // Compute partials.

//        Eigen::Vector3d partialsThirdOrderCoefficient =
//                2.0 * ( Eigen::Vector3d( ) <<
//                  instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 ),

//                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 ),

//                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 ) ).finished( ).transpose( );


//        Eigen::Vector3d partialsSecondOrderCoefficient =
//                3.0 * ( Eigen::Vector3d( ) <<
//                        instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 ),

//                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 ),

//                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 ) ).finished( ).transpose( );


//        Eigen::Vector3d partialsFirstOrderCoefficient =
//                2.0 * ( Eigen::Vector3d( ) <<
//                        instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 0 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 ),

//                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 1 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 ),

//                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 2 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
//                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 ) ).finished( ).transpose( );


//        Eigen::Vector3d partialsZeroOrderCoefficient =
//                2.0 * ( Eigen::Vector3d( ) <<
//                        instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 0 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 0 ),

//                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 1 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 1 ),

//                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 2 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 2 ) ).finished( ).transpose( );



//        // Assign partials to the right link leg.

//        if ( currentLegIndex == 0 ) // If first transmitter - receiver leg
//        {
//            partials.first.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient;
//            partials.first.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient;
//            partials.first.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient;
//            partials.first.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient;
//        }
//        else // If second transmitter - receiver leg
//        {
//            partials.second.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient;
//            partials.second.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient;
//            partials.second.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient;
//            partials.second.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient;
//        }


//    }

//    return partials;
//}


Eigen::Vector4d MutualApproximationScaling::computeCubicPolynomialCoefficients( )
{
    return ( Eigen::Vector4d( ) <<

              instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ],

            3.0 * ( instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),

            2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ]
            + instrumentalFrameRelativeVelocity_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativeVelocity_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),

            2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ) ).finished( );
}


void MutualApproximationScaling::computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition(
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtReceiverPosition )
{

    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition;
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition;
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition;

    // Successively compute partials w.r.t. position of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        if ( currentLinkEndIndex == 0 ) // Partials w.r.t. first transmitter position
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_;
            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_;
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // Partials w.r.t. second transmitter position
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_;
            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // Partials w.r.t. receiver position
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_;

            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_; /*- partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_.first
                    - partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_.second;*/

            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_; /*- partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_.first
                    - partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_.second;*/
        }


        // Compute partials.

        Eigen::Vector3d partialsThirdOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                  instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsSecondOrderCoefficient =
                3.0 * ( Eigen::Vector3d( ) <<
                        instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 ),

                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 ),

                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsFirstOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                        instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 0 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 0 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 1 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 1 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 2 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition( 1, 2 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsZeroOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                        instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 0 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 0 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 0 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 1 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 1 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 1 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 0, 2 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition( 1, 2 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition( 1, 2 ) ).finished( ); //.transpose( );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialOfCubicPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverPosition.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverPosition.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverPosition.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
    }
}


Eigen::Vector3d MutualApproximationScaling::computeDepressedCubicPolynomialCoefficients(  )
{
    return ( Eigen::Vector3d( ) <<
             cubicPolynomialCoefficients_[ 1 ] / cubicPolynomialCoefficients_[ 0 ],
            cubicPolynomialCoefficients_[ 2 ] / cubicPolynomialCoefficients_[ 0 ],
            cubicPolynomialCoefficients_[ 3 ] / cubicPolynomialCoefficients_[ 0 ] ).finished( );
}


void MutualApproximationScaling::computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianPosition(
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtReceiverPosition )
{
//    // Set partial matrix.
//    std::pair< Eigen::Matrix< double, 3, 3 >, Eigen::Matrix< double, 3, 3 > > partials
//            = std::make_pair( Eigen::Matrix< double, 3, 3 >::Zero( ), Eigen::Matrix< double, 3, 3 >::Zero( ) );


    // Compute coefficients of the cubic polynomial expression of the central instant t0.

//    Eigen::Vector4d cubicPolynomialCoefficients =
//            ( Eigen::Vector4d( ) <<

//              instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
//            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ],

//            3.0 * ( instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),

//            2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
//            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ]
//            + instrumentalFrameRelativeVelocity_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//            + instrumentalFrameRelativeVelocity_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),

//            2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ) ).finished( );


//    // Compute partials of cubic polynomial coefficients wrt link end position.
//    std::pair< Eigen::Matrix< double, 4, 3 >, Eigen::Matrix< double, 4, 3 > > partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition
//            = computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition( );

    /// NEW WAY TO COMPUTE PARTIALS OF CUBIC POLYNOMIAL COEFFICIENTS W.R.T. LINK END POSITION
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtReceiverPosition;
    computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition( partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
                                                                     partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
                                                                     partialOfCubicPolynomialCoefficientsWrtReceiverPosition );

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndPosition = partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndPosition = partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndPosition = partialOfCubicPolynomialCoefficientsWrtReceiverPosition;
        }

        Eigen::Matrix< double, 3, 3 > partials = Eigen::Matrix< double, 3, 3 >::Zero( );

        // Compute partials of the second order coefficient of the depressed cubic polynomial.
        partials.block( 0, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 1, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 1 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 ) );

//    partials.second.block( 0, 0, 1, 3 ) =
//            ( 1.0 / ( cubicPolynomialCoefficients[ 0 ] * cubicPolynomialCoefficients[ 0 ] ) )
//            * ( cubicPolynomialCoefficients[ 0 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 1, 0, 1, 3 )
//            - cubicPolynomialCoefficients[ 1 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) );



        // Compute partials of the first order coefficient of the depressed cubic polynomial.
        partials.block( 1, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 2, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 2 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 ) );

//        partials.second.block( 1, 0, 1, 3 ) =
//                ( 1.0 / ( cubicPolynomialCoefficients[ 0 ] * cubicPolynomialCoefficients[ 0 ] ) )
//                * ( cubicPolynomialCoefficients[ 0 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 2, 0, 1, 3 )
//                - cubicPolynomialCoefficients[ 2 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) );



        // Compute partials of the zero order coefficient of the depressed cubic polynomial.
        partials.block( 2, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 3, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 3 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 ) );

//        partials.second.block( 2, 0, 1, 3 ) =
//                ( 1.0 / ( cubicPolynomialCoefficients[ 0 ] * cubicPolynomialCoefficients[ 0 ] ) )
//                * ( cubicPolynomialCoefficients[ 0 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 3, 0, 1, 3 )
//                - cubicPolynomialCoefficients[ 3 ] * partialsOfCubicPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 ) );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterPosition = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtReceiverPosition = partials;
        }

    }



//    return partials;
}


//double computePartialYwrtTime( Eigen::Vector3d cartesianPositionFirstObject,
//                               Eigen::Vector3d cartesianPositionSecondObject,
//                               Eigen::Vector3d cartesianVelocityFirstObject,
//                               Eigen::Vector3d cartesianVelocitySecondObject )
//{
//    double partialDeclinationFirstObject = computePartialDeclinationWrtTime( cartesianPositionFirstObject,
//                                                                             cartesianVelocityFirstObject );
//    double partialDeclinationSecondObject = computePartialDeclinationWrtTime( cartesianPositionSecondObject,
//                                                                              cartesianVelocitySecondObject );

//    double partial = partialDeclinationSecondObject - partialDeclinationFirstObject;

//    return partial;
//}


void MutualApproximationScaling::computePartialsOfCentralInstantWrtLinkEndPosition(
        const std::vector< double > times )
{

    double intermediateVariableQ = ( 3.0 * depressedCubicPolynomialCoefficients_[ 1 ]
            - depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] ) / 9.0;
    double intermediateVariableR = ( 9.0 * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 1 ]
            - 27.0 * depressedCubicPolynomialCoefficients_[ 2 ]
            - 2.0 * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] ) / 54.0;

    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition;
    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition;
    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverPosition;
    computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianPosition( partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
                                                                              partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition,
                                                                              partialOfDepressedPolynomialCoefficientsWrtReceiverPosition );

    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtFirstTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtSecondTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtReceiverPosition;
    computePartialOfIntermediateVariableQWrtLinkEndPosition(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
                partialsIntermediateVariableQwrtFirstTransmitterPosition, partialsIntermediateVariableQwrtSecondTransmitterPosition,
                partialsIntermediateVariableQwrtReceiverPosition );

    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtFirstTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtSecondTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtReceiverPosition;
    computePartialOfIntermediateVariableRWrtLinkEndPosition(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
                partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
                partialsIntermediateVariableRwrtReceiverPosition );

    if ( intermediateVariableBeta < 0 )
    {
        std::cout << "BETA NEGATITVE" << "\n\n";

        double thetaAngle = computeAngleThetaRealSolutionsCubicEquation( intermediateVariableQ, intermediateVariableR );

        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtReceiverPosition;
        computePartialOfAngleThetaCubicEquationWrtLinkEndPosition( intermediateVariableQ, intermediateVariableR,
                                                                   partialsIntermediateVariableQwrtFirstTransmitterPosition,
                                                                   partialsIntermediateVariableQwrtSecondTransmitterPosition,
                                                                   partialsIntermediateVariableQwrtReceiverPosition,
                                                                   partialsIntermediateVariableRwrtFirstTransmitterPosition,
                                                                   partialsIntermediateVariableRwrtSecondTransmitterPosition,
                                                                   partialsIntermediateVariableRwrtReceiverPosition,
                                                                   partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition,
                                                                   partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition,
                                                                   partialOfAngleThetaCubicEquationWrtReceiverPosition );

        double firstSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( thetaAngle / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;
        double secondSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;
        double thirdSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;

        std::cout << "estimated central instant: " << times[ 0 ] << "\n\n";
        std::cout << "firstSolutionCentralInstant: " << firstSolutionCentralInstant << "\n\n";
        std::cout << "secondSolutionCentralInstant: " << secondSolutionCentralInstant << "\n\n";
        std::cout << "thirdSolutionCentralInstant: " << thirdSolutionCentralInstant << "\n\n";


        // Compute partials of central instant w.r.t. link ends position.

        if ( ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) &&
             ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
        {
            std::cout << "FIRST SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverPosition_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
        }
        else if ( ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
                  ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
        {
            std::cout << "SECOND SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverPosition_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
        }
        else if ( ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
                  ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) )
        {
            std::cout << "THIRD SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverPosition_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
        }

    }

    else
    {

        std::cout << "BETA POSITIVE" << "\n\n";

        double intermediateVariableS = 0.0;
        if ( ( intermediateVariableR + sqrt( intermediateVariableBeta ) ) >= 0 )
        {
            intermediateVariableS = std::pow( intermediateVariableR + sqrt( intermediateVariableBeta ), 1.0 / 3.0 );
        }
        else
        {
            intermediateVariableS = - std::pow( std::fabs( intermediateVariableR + std::sqrt( intermediateVariableBeta ) ), 1.0 / 3.0 );
        }

        double intermediateVariableT = 0.0;
        if ( ( intermediateVariableR - sqrt( intermediateVariableBeta ) ) >= 0 )
        {
            intermediateVariableT = std::pow( intermediateVariableR - sqrt( intermediateVariableBeta ), 1.0 / 3.0 );
        }
        else
        {
            intermediateVariableT = - std::pow( std::fabs( intermediateVariableR - std::sqrt( intermediateVariableBeta ) ), 1.0 / 3.0 );
        }

        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtFirstTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtSecondTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtReceiverPosition;
        computePartialOfIntermediateVariableTWrtLinkEndPosition(
                    intermediateVariableQ, intermediateVariableR, intermediateVariableT, partialsIntermediateVariableQwrtFirstTransmitterPosition,
                    partialsIntermediateVariableQwrtSecondTransmitterPosition, partialsIntermediateVariableQwrtReceiverPosition,
                    partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
                    partialsIntermediateVariableRwrtReceiverPosition, partialsIntermediateVariableTwrtFirstTransmitterPosition,
                    partialsIntermediateVariableTwrtSecondTransmitterPosition, partialsIntermediateVariableTwrtReceiverPosition );

        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtFirstTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtSecondTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtReceiverPosition;
        computePartialOfIntermediateVariableSWrtLinkEndPosition(
                    intermediateVariableQ, intermediateVariableR, intermediateVariableS, partialsIntermediateVariableQwrtFirstTransmitterPosition,
                    partialsIntermediateVariableQwrtSecondTransmitterPosition, partialsIntermediateVariableQwrtReceiverPosition,
                    partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
                    partialsIntermediateVariableRwrtReceiverPosition, partialsIntermediateVariableSwrtFirstTransmitterPosition,
                    partialsIntermediateVariableSwrtSecondTransmitterPosition, partialsIntermediateVariableSwrtReceiverPosition );

        // Compute partials of central instant t0 w.r.t. link end positon.
        partialsOfCentralInstantWrtFirstTransmitterPosition_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtFirstTransmitterPosition + partialsIntermediateVariableTwrtFirstTransmitterPosition;

        partialsOfCentralInstantWrtSecondTransmitterPosition_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtSecondTransmitterPosition + partialsIntermediateVariableTwrtSecondTransmitterPosition;

        partialsOfCentralInstantWrtReceiverPosition_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtReceiverPosition + partialsIntermediateVariableTwrtReceiverPosition;
    }
}


//! Function to check whether all the required dependent variables are included in the interface.
void MutualApproximationScaling::checkRequiredDependentVariablesFromInterface( const observation_models::LinkEnds linkEnds )
{
    if ( dependentVariablesInterfaceIdsAndIndices_.count(  "Total acceleration in inertial frame of " + linkEnds.at( observation_models::receiver ).first ) == 0 )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration of receiver in the dependent variables interface." );
    }
    if ( ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                           + linkEnds.at( observation_models::receiver ).first
                                                           + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                              + linkEnds.at( observation_models::receiver ).first
                                                              + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                                 + linkEnds.at( observation_models::receiver ).first
                                                                 + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first ) == 0 ) )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration partial of receiver w.r.t. link ends translational states"
                                  " in the dependent variables interface." );
    }
    if ( dependentVariablesInterfaceIdsAndIndices_.count(  "Total acceleration in inertial frame of " + linkEnds.at( observation_models::transmitter ).first ) == 0 )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration of transmitter in the dependent variables interface." );
    }
    if ( ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                           + linkEnds.at( observation_models::transmitter ).first
                                                           + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                              + linkEnds.at( observation_models::transmitter ).first
                                                              + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                              + linkEnds.at( observation_models::transmitter ).first
                                                              + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first ) == 0 ) )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration partial of transmitter w.r.t. link ends translational states"
                                  " in the dependent variables interface." );
    }
    if ( dependentVariablesInterfaceIdsAndIndices_.count(  "Total acceleration in inertial frame of " + linkEnds.at( observation_models::transmitter2 ).first ) == 0 )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration of transmitter2 in the dependent variables interface." );
    }
    if ( ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                           + linkEnds.at( observation_models::transmitter2 ).first
                                                           + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                              + linkEnds.at( observation_models::transmitter2 ).first
                                                              + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first ) == 0 )
         || ( dependentVariablesInterfaceIdsAndIndices_.count( "Total acceleration partial w.r.t. translational state, for acceleration of "
                                                              + linkEnds.at( observation_models::transmitter2 ).first
                                                              + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first ) == 0 ) )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find total acceleration partial of transmitter2 w.r.t. link ends translational states"
                                  " in the dependent variables interface." );
    }
}

//! Function to retrieve the relative accelerations of the two transmitters, and the associated acceleration partials
//!  from the dependent variables interface.
void MutualApproximationScaling::retrieveRelativeAccelerationsAndAssociatedPartialsFromDependentVariables(
        const std::vector< double >& times, const observation_models::LinkEnds linkEnds )
{
    // Compute relative accelerations of the two transmitters w.r.t the common receiver.
    cartesianAccelerationFirstTransmitterWrtReceiver_ = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration in inertial frame of " + linkEnds.at( observation_models::transmitter ).first, 3, times[ 0 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration in inertial frame of " + linkEnds.at( observation_models::receiver ).first, 3, times[ 2 ] );

    cartesianAccelerationSecondTransmitterWrtReceiver_ = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration in inertial frame of " + linkEnds.at( observation_models::transmitter2 ).first, 3, times[ 1 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration in inertial frame of " + linkEnds.at( observation_models::receiver ).first, 3, times[ 2 ] );

    // Compute partials of the relative accelerations of the two transmitters w.r.t. the common receiver.
    Eigen::VectorXd partialAccelerationFirstTransmitterWrtReceiverPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first, 18, times[ 0 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationFirstTransmitterWrtReceiverPosition_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtReceiverPosition.segment( i * 6, 3 ).transpose( );
    }

    Eigen::VectorXd partialAccelerationFirstTransmitterWrtTransmitterPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first, 18, times[ 0 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationFirstTransmitterWrtTransmitterPosition_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtTransmitterPosition.segment( i * 6, 3 ).transpose( );
    }

    Eigen::VectorXd partialAccelerationFirstTransmitterWrtOtherTransmitterPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first, 18, times[ 0 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtOtherTransmitterPosition.segment( i * 6, 3 ).transpose( );
    }

    Eigen::VectorXd partialAccelerationSecondTransmitterWrtReceiverPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter2 ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first, 18, times[ 1 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::receiver ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationSecondTransmitterWrtReceiverPosition_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtReceiverPosition.segment( i * 6, 3 ).transpose( );
    }

    Eigen::VectorXd partialAccelerationSecondTransmitterWrtTransmitterPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter2 ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first, 18, times[ 1 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter2 ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationSecondTransmitterWrtTransmitterPosition_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtTransmitterPosition.segment( i * 6, 3 ).transpose( );
    }

    Eigen::VectorXd partialAccelerationSecondTransmitterWrtOtherTransmitterPosition = dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::transmitter2 ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first, 18, times[ 0 ] )
            - dependentVariablesInterface_->getSingleDependentVariable(
                "Total acceleration partial w.r.t. translational state, for acceleration of " + linkEnds.at( observation_models::receiver ).first
                + " w.r.t. state of " + linkEnds.at( observation_models::transmitter ).first, 18, times[ 2 ] );
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtOtherTransmitterPosition.segment( i * 6, 3 ).transpose( );
    }
}

//! Update the scaling object to the current times and states
void MutualApproximationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                         const std::vector< double >& times,
                                         const observation_models::LinkEndType fixedLinkEnd,
                                         const observation_models::LinkEnds linkEnds,
                                         const Eigen::VectorXd currentObservation )
{
    std::cout << "BEGINNING UPDATE FUNCTION IN MUTUAL APPROXIMATION SCALING" << "\n\n";

    if ( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating the mutual approximation scaling object, "
                                  "fixed link end time different from receiver." );
    }

    if( ( linkEnds.count( observation_models::receiver ) == 0 ) ||
            ( linkEnds.count( observation_models::transmitter ) == 0 ) ||
            ( linkEnds.count( observation_models::transmitter2 ) == 0 ) )
    {
        throw std::runtime_error( "Error when updating mutual approximation partials scaling, did not find transmitter, transmitter2 or receiver in link ends" );
    }

    // Check that the required dependent variables are provided by the dependent variables interface.
    checkRequiredDependentVariablesFromInterface( linkEnds );

    // Retrieve the relative accelerations of the two transmitters, and the associated acceleration partials.
    retrieveRelativeAccelerationsAndAssociatedPartialsFromDependentVariables( times, linkEnds );

    // Retrieve states of the two transmitters and the receiver.
    firstTransmitterState_ = linkEndStates[ 0 ];
    secondTransmitterState_ = linkEndStates[ 1 ];
    receiverState_ = linkEndStates[ 2 ];

    Eigen::Vector3d relativeRangeVectorFirstTransmitter = ( receiverState_ - firstTransmitterState_ ).segment( 0, 3 );
    Eigen::Vector3d relativeRangeVectorSecondTransmitter = ( receiverState_ - secondTransmitterState_ ).segment( 0, 3 );

//    std::cout << "relativeRangeVectorFirstTransmitter: " << relativeRangeVectorFirstTransmitter.transpose( ) << "\n\n";
//    std::cout << "relativeRangeVectorSecondTransmitter: " << relativeRangeVectorSecondTransmitter.transpose( ) << "\n\n";

    // Compute right ascension and declination of the first and second transmitters as seen from the receiver.
    std::pair< double, double > rightAscensionAndDeclinationFirstTransmitter =
            computeRightAscensionAndDeclination( relativeRangeVectorFirstTransmitter );
    std::pair< double, double > rightAscensionAndDeclinationSecondTransmitter =
            computeRightAscensionAndDeclination( relativeRangeVectorSecondTransmitter );

    rightAscensionFirstTransmitter_ = rightAscensionAndDeclinationFirstTransmitter.first;
    declinationFirstTransmitter_ = rightAscensionAndDeclinationFirstTransmitter.second;
    rightAscensionSecondTransmitter_ = rightAscensionAndDeclinationSecondTransmitter.first;
    declinationSecondTransmitter_ = rightAscensionAndDeclinationSecondTransmitter.second;

    averageDeclination_ = ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0;


    // Compute first partials of right ascension and declination of the two transmitters wrt time.
    partialOfRightAscensionFirstTransmitterWrtTime_ = computePartialOfRightAscensionWrtTime(
                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

    partialOfRightAscensionSecondTransmitterWrtTime_ = computePartialOfRightAscensionWrtTime(
                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

    partialOfDeclinationFirstTransmitterWrtTime_ = computePartialOfDeclinationWrtTime(
                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

    partialOfDeclinationSecondTransmitterWrtTime_ = computePartialOfDeclinationWrtTime(
                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ) );

//    std::cout << "first partial right ascension w.r.t. time: " << partialOfRightAscensionSecondTransmitterWrtTime_ << "\n\n";
//    std::cout << "first partial declination w.r.t. time: " << partialOfDeclinationSecondTransmitterWrtTime_ << "\n\n";


    // Compute second partials of right ascension and declination of the two transmitters wrt time.
    secondPartialOfRightAscensionFirstTransmitterWrtTime_ = computeSecondPartialRightAscensionWrtTime(
                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
                cartesianAccelerationFirstTransmitterWrtReceiver_ );

    secondPartialOfRightAscensionSecondTransmitterWrtTime_ = computeSecondPartialRightAscensionWrtTime(
                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ),
                cartesianAccelerationSecondTransmitterWrtReceiver_ );

    secondPartialOfDeclinationFirstTransmitterWrtTime_ = computeSecondPartialDeclinationWrtTime(
                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
                cartesianAccelerationFirstTransmitterWrtReceiver_ );

    secondPartialOfDeclinationSecondTransmitterWrtTime_ = computeSecondPartialDeclinationWrtTime(
                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ),
                cartesianAccelerationSecondTransmitterWrtReceiver_ );


    // Compute relative position, velocity and acceleration between two transmitters in the instrumental frame of
    // the receiver.
    instrumentalFrameRelativePosition_ = computeRelativePositionInInstrumentalFrame( );

//    std::cout << "instrumental frame relative position: " << instrumentalFrameRelativePosition_.transpose( ) << "\n\n";

    instrumentalFrameRelativeVelocity_ = computeRelativeVelocityInInstrumentalFrame( );

//    std::cout << "instrumental frame relative velocity: " << instrumentalFrameRelativeVelocity_.transpose( ) << "\n\n";

    instrumentalFrameRelativeAcceleration_ = computeRelativeAccelerationInInstrumentalFrame(
                cartesianAccelerationFirstTransmitterWrtReceiver_,
                cartesianAccelerationSecondTransmitterWrtReceiver_ );

//    std::cout << "instrumental frame relative acceleration: " << instrumentalFrameRelativeAcceleration_.transpose( ) << "\n\n";



    // Compute the coefficients of the depressed cubic polynomial expressing the central instant t0 as a function
    // of the relative position, velocity and acceleration between the two transmitters in the instrumental frame
    // of the receiver, at t0.

//    double squaredNormInstrumentalFrameRelativeAcceleration =
//            instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
//            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ];

//    Eigen::Vector3d coefficientsDepressedCubicPolynomialCentralInstant = ( Eigen::Vector3d( ) <<
//              3.0 * ( instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//              + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),
//              2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
//            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ]
//            + instrumentalFrameRelativeVelocity_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//            + instrumentalFrameRelativeVelocity_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ),
//              2.0 * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
//            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] ) ).finished( );

//    coefficientsDepressedCubicPolynomialCentralInstant /=
//            ( instrumentalFrameRelativeAcceleration_[ 0 ] * instrumentalFrameRelativeAcceleration_[ 0 ]
//            + instrumentalFrameRelativeAcceleration_[ 1 ] * instrumentalFrameRelativeAcceleration_[ 1 ] );



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian position.
    /*partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_ =*/ computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPosition( );
    /*partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_ =*/ computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPosition( );

    // Compute partials of relative position between the two transmitters in the instrumental frame wrt link end cartesian position.


    // Compute partials of relative velocity between the two transmitters in the instrumental frame wrt link end cartesian position.


//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition: " << partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition_.first.transpose( ) << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition: " << partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition_.first.transpose( ) << "\n\n";



//    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian positions.
//    partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndPosition_ = computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
//                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_,
//                partialAccelerationFirstTransmitterWrtLinkEndPosition_, partialAccelerationSecondTransmitterWrtLinkEndPosition_ );

    // partials of relative acceleration between two transmitters in the instrumental frame w.r.t. link end positions (returned by reference).
    /// DONE DIFFERENTLY THAN THE ABOVE
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_,
                partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_, partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_,
                partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ );

//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ << "\n\n";


//    std::cout << "----------------------------------------------------------------------------------------------------------------------" << "\n\n";
//    std::cout << "--------------- TEST POSITION PARTIALS OF INTERMEDIATE VARIABLES USED TO DERIVE THE CENTRAL INSTANT ------------------" << "\n\n";
//    std::cout << "----------------------------------------------------------------------------------------------------------------------" << "\n\n";

    // Compute coefficients of the cubic polynomial for the central instant t0.
    cubicPolynomialCoefficients_ = computeCubicPolynomialCoefficients( );

    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
    depressedCubicPolynomialCoefficients_ = computeDepressedCubicPolynomialCoefficients( );

    // Compute partials of central instant w.r.t. link ends positions.
    computePartialsOfCentralInstantWrtLinkEndPosition( times );



//    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
//    Eigen::Vector3d depressedCubicPolynomialCoefficients =
//            ( Eigen::Vector3d( ) <<
//              cubicPolynomialCoefficients[ 1 ] / cubicPolynomialCoefficients[ 0 ],
//            cubicPolynomialCoefficients[ 2 ] / cubicPolynomialCoefficients[ 0 ],
//            cubicPolynomialCoefficients[ 3 ] / cubicPolynomialCoefficients[ 0 ] ).finished( );

//    double intermediateVariableQ = ( 3.0 * depressedCubicPolynomialCoefficients[ 1 ]
//            - depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ] ) / 9.0;
//    double intermediateVariableR = ( 9.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 1 ]
//            - 27.0 * depressedCubicPolynomialCoefficients[ 2 ]
//            - 2.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ] ) / 54.0;

////    std::cout << "intermediate variable Q (model): " << intermediateVariableQ << "\n\n";
////    std::cout << "intermediate variable R (model): " << intermediateVariableR << "\n\n";

//    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
//            + intermediateVariableR * intermediateVariableR;

//    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition;
//    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition;
//    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverPosition;
//    computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianPosition( partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
//                                                                              partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition,
//                                                                              partialOfDepressedPolynomialCoefficientsWrtReceiverPosition );

////    std::cout << "partials depressed polynomial coefficients w.r.t. first transmitter position (model): " << partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition << "\n\n";
////    std::cout << "partials depressed polynomial coefficients w.r.t. second transmitter position (model): " << partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition << "\n\n";
////    std::cout << "partials depressed polynomial coefficients w.r.t. receiver position (model): " << partialOfDepressedPolynomialCoefficientsWrtReceiverPosition << "\n\n";

////    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableQwrtLinkEndPosition
////            = computePartialOfIntermediateVariableQWrtLinkEndPosition( depressedCubicPolynomialCoefficients,
////                                                                       partialDepressedPolynomialCoefficientsWrtLinkEndPosition );

//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtFirstTransmitterPosition;
//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtSecondTransmitterPosition;
//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtReceiverPosition;
//    computePartialOfIntermediateVariableQWrtLinkEndPosition(
//                depressedCubicPolynomialCoefficients, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
//                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
//                partialsIntermediateVariableQwrtFirstTransmitterPosition, partialsIntermediateVariableQwrtSecondTransmitterPosition,
//                partialsIntermediateVariableQwrtReceiverPosition );

////    std::cout << "partials intermediate Q w.r.t. first transmitter position (model): " << partialsIntermediateVariableQwrtFirstTransmitterPosition << "\n\n";
////    std::cout << "partials intermediate Q w.r.t. second transmitter position (model): " << partialsIntermediateVariableQwrtSecondTransmitterPosition << "\n\n";
////    std::cout << "partials intermediate Q w.r.t. receiver position (model): " << partialsIntermediateVariableQwrtReceiverPosition << "\n\n";

////    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableRwrtLinkEndPosition
////            = computePartialOfIntermediateVariableRWrtLinkEndPosition( depressedCubicPolynomialCoefficients,
////                                                                       partialDepressedPolynomialCoefficientsWrtLinkEndPosition );

//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtFirstTransmitterPosition;
//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtSecondTransmitterPosition;
//    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtReceiverPosition;
//    computePartialOfIntermediateVariableRWrtLinkEndPosition(
//                depressedCubicPolynomialCoefficients, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
//                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
//                partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
//                partialsIntermediateVariableRwrtReceiverPosition );

////    std::cout << "partials intermediate R w.r.t. first transmitter position (model): " << partialsIntermediateVariableRwrtFirstTransmitterPosition << "\n\n";
////    std::cout << "partials intermediate R w.r.t. second transmitter position (model): " << partialsIntermediateVariableRwrtSecondTransmitterPosition << "\n\n";
////    std::cout << "partials intermediate R w.r.t. receiver position (model): " << partialsIntermediateVariableRwrtReceiverPosition << "\n\n";


////    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsOfCentralInstantWrtLinkEndPosition;

////    std::cout << "intermediate variable beta: " << intermediateVariableBeta << "\n\n";

//    if ( intermediateVariableBeta < 0 )
//    {
//        std::cout << "BETA NEGATITVE" << "\n\n";

//        double thetaAngle = computeAngleThetaRealSolutionsCubicEquation( intermediateVariableQ, intermediateVariableR );

////        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsOfAngleThetaCubicEquationWrtLinkEndPosition
////                = computePartialOfAngleThetaCubicEquationWrtLinkEndPosition( intermediateVariableQ, intermediateVariableR,
////                                                                             partialsIntermediateVariableQwrtLinkEndPosition,
////                                                                             partialsIntermediateVariableRwrtLinkEndPosition );

//        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtReceiverPosition;
//        computePartialOfAngleThetaCubicEquationWrtLinkEndPosition( intermediateVariableQ, intermediateVariableR,
//                                                                   partialsIntermediateVariableQwrtFirstTransmitterPosition,
//                                                                   partialsIntermediateVariableQwrtSecondTransmitterPosition,
//                                                                   partialsIntermediateVariableQwrtReceiverPosition,
//                                                                   partialsIntermediateVariableRwrtFirstTransmitterPosition,
//                                                                   partialsIntermediateVariableRwrtSecondTransmitterPosition,
//                                                                   partialsIntermediateVariableRwrtReceiverPosition,
//                                                                   partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition,
//                                                                   partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition,
//                                                                   partialOfAngleThetaCubicEquationWrtReceiverPosition );

////        std::cout << "partials angle theta cubic equation w.r.t. first transmitter position (model): " << partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition << "\n\n";
////        std::cout << "partials angle theta cubic equation w.r.t. second transmitter position (model): " << partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition << "\n\n";
////        std::cout << "partials angle theta cubic equation w.r.t. receiver position (model): " << partialOfAngleThetaCubicEquationWrtReceiverPosition << "\n\n";


//        double firstSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( thetaAngle / 3.0 )
//                - depressedCubicPolynomialCoefficients[ 0 ] / 3.0;
//        double secondSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
//                - depressedCubicPolynomialCoefficients[ 0 ] / 3.0;
//        double thirdSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
//                - depressedCubicPolynomialCoefficients[ 0 ] / 3.0;

//        std::cout << "estimated central instant: " << times[ 0 ] << "\n\n";
//        std::cout << "firstSolutionCentralInstant: " << firstSolutionCentralInstant << "\n\n";
//        std::cout << "secondSolutionCentralInstant: " << secondSolutionCentralInstant << "\n\n";
//        std::cout << "thirdSolutionCentralInstant: " << thirdSolutionCentralInstant << "\n\n";


//        // Compute partials of central instant w.r.t. link end states

////        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtFirstTransmitterPosition;
////        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtSecondTransmitterPosition;
////        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtReceiverPosition;


////        partialsOfCentralInstantWrtLinkEndPosition.first = std::cos( thetaAngle / 3.0 ) / sqrt( - intermediateVariableQ )
////                * ( - partialsIntermediateVariableQwrtLinkEndPosition.first )
////                - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( thetaAngle / 3.0 ) * partialsOfAngleThetaCubicEquationWrtLinkEndPosition.first
////                - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 0, 0, 1, 3 );

////        partialsOfCentralInstantWrtLinkEndPosition.second = std::cos( thetaAngle / 3.0 ) / sqrt( - intermediateVariableQ )
////                * ( - partialsIntermediateVariableQwrtLinkEndPosition.second )
////                - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( thetaAngle / 3.0 ) * partialsOfAngleThetaCubicEquationWrtLinkEndPosition.second
////                - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 );

//        if ( ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) &&
//             ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
//        {
//            std::cout << "FIRST SOLUTION CENTRAL INSTANT" << "\n\n";
//            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
//                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
//                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtReceiverPosition_ =
//                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
//        }
//        else if ( ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
//                  ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
//        {
//            std::cout << "SECOND SOLUTION CENTRAL INSTANT" << "\n\n";
//            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
//                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
//                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtReceiverPosition_ =
//                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
//        }
//        else if ( ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
//                  ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) )
//        {
//            std::cout << "THIRD SOLUTION CENTRAL INSTANT" << "\n\n";
//            partialsOfCentralInstantWrtFirstTransmitterPosition_ =
//                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtFirstTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtSecondTransmitterPosition_ =
//                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtSecondTransmitterPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 );

//            partialsOfCentralInstantWrtReceiverPosition_ =
//                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
//                    * ( - partialsIntermediateVariableQwrtReceiverPosition )
//                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
//                    * partialOfAngleThetaCubicEquationWrtReceiverPosition
//                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 );
//        }

////        std::cout << "partials central instant w.r.t. first transmitter position (model): " << partialsOfCentralInstantWrtFirstTransmitterPosition_ << "\n\n";
////        std::cout << "partials central instant w.r.t. second transmitter position (model): " << partialsOfCentralInstantWrtSecondTransmitterPosition_ << "\n\n";
////        std::cout << "partials central instant w.r.t. receiver position (model): " << partialsOfCentralInstantWrtReceiverPosition_ << "\n\n";

////        partialsOfCentralInstantWrtLinkEndPosition.second = std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
////                * ( - partialsIntermediateVariableQwrtLinkEndPosition.second )
////                - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI) / 3.0 )
////                * partialsOfAngleThetaCubicEquationWrtLinkEndPosition.second
////                - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 );


////        partialsOfCentralInstantWrtLinkEndPosition.first = std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
////                * ( - partialsIntermediateVariableQwrtLinkEndPosition.first )
////                - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
////                * partialsOfAngleThetaCubicEquationWrtLinkEndPosition.first
////                - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.first.block( 0, 0, 1, 3 );

////        partialsOfCentralInstantWrtLinkEndPosition.second = std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
////                * ( - partialsIntermediateVariableQwrtLinkEndPosition.second )
////                - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
////                * partialsOfAngleThetaCubicEquationWrtLinkEndPosition.second
////                - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 );
//    }

//    else
//    {

//        std::cout << "BETA POSITIVE" << "\n\n";

//        double intermediateVariableS = 0.0;
//        if ( ( intermediateVariableR + sqrt( intermediateVariableBeta ) ) >= 0 )
//        {
//            intermediateVariableS = std::pow( intermediateVariableR + sqrt( intermediateVariableBeta ), 1.0 / 3.0 );
//        }
//        else
//        {
//            intermediateVariableS = - std::pow( std::fabs( intermediateVariableR + std::sqrt( intermediateVariableBeta ) ), 1.0 / 3.0 );
//        }

//        double intermediateVariableT = 0.0;
//        if ( ( intermediateVariableR - sqrt( intermediateVariableBeta ) ) >= 0 )
//        {
//            intermediateVariableT = std::pow( intermediateVariableR - sqrt( intermediateVariableBeta ), 1.0 / 3.0 );
//        }
//        else
//        {
//            intermediateVariableT = - std::pow( std::fabs( intermediateVariableR - std::sqrt( intermediateVariableBeta ) ), 1.0 / 3.0 );
//        }

////        std::cout << "intermediate variable S (model): " << intermediateVariableS << "\n\n";
////        std::cout << "intermediate variable T (model): " << intermediateVariableT << "\n\n";


//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtFirstTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtSecondTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtReceiverPosition;
//        computePartialOfIntermediateVariableTWrtLinkEndPosition(
//                    intermediateVariableQ, intermediateVariableR, intermediateVariableT, partialsIntermediateVariableQwrtFirstTransmitterPosition,
//                    partialsIntermediateVariableQwrtSecondTransmitterPosition, partialsIntermediateVariableQwrtReceiverPosition,
//                    partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
//                    partialsIntermediateVariableRwrtReceiverPosition, partialsIntermediateVariableTwrtFirstTransmitterPosition,
//                    partialsIntermediateVariableTwrtSecondTransmitterPosition, partialsIntermediateVariableTwrtReceiverPosition );

////        std::cout << "partials T w.r.t. first transmitter position (model): " << partialsIntermediateVariableTwrtFirstTransmitterPosition << "\n\n";
////        std::cout << "partials T w.r.t. second transmitter position (model): " << partialsIntermediateVariableTwrtSecondTransmitterPosition << "\n\n";
////        std::cout << "partials T w.r.t. receiver (model): " << partialsIntermediateVariableTwrtReceiverPosition << "\n\n";

////        std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > partialsIntermediateVariableSwrtLinkEndPosition =
////                computePartialOfIntermediateVariableSWrtLinkEndPosition( intermediateVariableQ, intermediateVariableR,
////                                                                         partialsIntermediateVariableQwrtLinkEndPosition,
////                                                                         partialsIntermediateVariableRwrtLinkEndPosition );

//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtFirstTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtSecondTransmitterPosition;
//        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtReceiverPosition;
//        computePartialOfIntermediateVariableSWrtLinkEndPosition(
//                    intermediateVariableQ, intermediateVariableR, intermediateVariableS, partialsIntermediateVariableQwrtFirstTransmitterPosition,
//                    partialsIntermediateVariableQwrtSecondTransmitterPosition, partialsIntermediateVariableQwrtReceiverPosition,
//                    partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
//                    partialsIntermediateVariableRwrtReceiverPosition, partialsIntermediateVariableSwrtFirstTransmitterPosition,
//                    partialsIntermediateVariableSwrtSecondTransmitterPosition, partialsIntermediateVariableSwrtReceiverPosition );

////        std::cout << "partials S w.r.t. first transmitter position (model): " << partialsIntermediateVariableSwrtFirstTransmitterPosition << "\n\n";
////        std::cout << "partials S w.r.t. second transmitter position (model): " << partialsIntermediateVariableSwrtSecondTransmitterPosition << "\n\n";
////        std::cout << "partials S w.r.t. receiver (model): " << partialsIntermediateVariableSwrtReceiverPosition << "\n\n";


//        // Compute partials of central instant t0 w.r.t. link end positon.
//        partialsOfCentralInstantWrtFirstTransmitterPosition_ =
//                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition.block( 0, 0, 1, 3 )
//                + partialsIntermediateVariableSwrtFirstTransmitterPosition + partialsIntermediateVariableTwrtFirstTransmitterPosition;

//        partialsOfCentralInstantWrtSecondTransmitterPosition_ =
//                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition.block( 0, 0, 1, 3 )
//                + partialsIntermediateVariableSwrtSecondTransmitterPosition + partialsIntermediateVariableTwrtSecondTransmitterPosition;

//        partialsOfCentralInstantWrtReceiverPosition_ =
//                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverPosition.block( 0, 0, 1, 3 )
//                + partialsIntermediateVariableSwrtReceiverPosition + partialsIntermediateVariableTwrtReceiverPosition;

////        std::cout << "partials central instant w.r.t. first transmitter position (model): " << partialsOfCentralInstantWrtFirstTransmitterPosition_ << "\n\n";
////        std::cout << "partials central instant w.r.t. second transmitter position (model): " << partialsOfCentralInstantWrtSecondTransmitterPosition_ << "\n\n";
////        std::cout << "partials central instant w.r.t. receiver position (model): " << partialsOfCentralInstantWrtReceiverPosition_ << "\n\n";

////        partialsOfCentralInstantWrtLinkEndPosition.second = - 1.0 / 3.0 * partialDepressedPolynomialCoefficientsWrtLinkEndPosition.second.block( 0, 0, 1, 3 )
////                + partialsIntermediateVariableSwrtLinkEndPosition.first
////                + partialsIntermediateVariableTwrtLinkEndPosition.first;

//        std::cout << "END UPDATE FUNCTION IN MUTUAL APPROXIMATION SCALING" << "\n\n";

//    }


////    // Compute reference scaling factors for angular position of both transmitters (at receiver).
////    std::vector< Eigen::Vector6d > linkEndStatesFirstTransmitterReceiver;
////    linkEndStatesFirstTransmitterReceiver.push_back( firstTransmitterState );
////    linkEndStatesFirstTransmitterReceiver.push_back( receiverState );

////    std::vector< Eigen::Vector6d > linkEndStatesSecondTransmitterReceiver;
////    linkEndStatesSecondTransmitterReceiver.push_back( secondTransmitterState );
////    linkEndStatesSecondTransmitterReceiver.push_back( receiverState );

////    referenceAngularPositionScaler_->update( linkEndStatesFirstTransmitterReceiver, times, fixedLinkEnd, currentObservation );
////    angularPositionScalingFactorFirstTransmitter_ = referenceAngularPositionScaler_->getScalingFactor( observation_models::receiver );
////    angularPositionLightTimeCorrectionScalingFirstTransmitter_ = referenceAngularPositionScaler_->getLightTimePartialScalingFactor( );

////    referenceAngularPositionScaler_->update( linkEndStatesSecondTransmitterReceiver, times, fixedLinkEnd, currentObservation );
////    angularPositionScalingFactorSecondTransmitter_ = referenceAngularPositionScaler_->getScalingFactor( observation_models::receiver );
////    angularPositionLightTimeCorrectionScalingSecondTransmitter_ = referenceAngularPositionScaler_->getLightTimePartialScalingFactor( );


////    //


//////    // Compute common scaling factor
//////    scalingFactor_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVectorFirstTransmitter, true );




//////    // Compute scaling for receiver reference
//////    if( fixedLinkEnd == observation_models::receiver )
//////    {
//////        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
//////                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
//////        referenceScalingFactor_ =
//////                scalingFactor_ *
//////                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
//////                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
//////    }
//////    // Compute scaling for transmitter reference
//////    else if( fixedLinkEnd == observation_models::transmitter )
//////    {
//////        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
//////                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
//////        referenceScalingFactor_ =
//////                scalingFactor_ *
//////                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
//////                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
//////    }


////    // TO BE MODIFIED, WORKS AS WELL FOR PARTIALS WRT LIGHT TIME CORRECTION PARAMETERS
////    scalingFactorXCoefficient_ = std::make_pair(
////                std::cos( averageDeclination_ ), - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
////                * std::sin( averageDeclination_ ) );
////    scalingFactorYCoefficient_ = 1.0;

////    XscalingFactors_ = //scalingFactorsXrightAscensionContribution_ =
////            std::make_pair( - angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) * scalingFactorXCoefficient_.first
////                            + angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorXCoefficient_.second,
////                            angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 ) * scalingFactorXCoefficient_.first
////                            + angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorXCoefficient_.second );

////    XlightTimeCorrectionScalingFactors_ = //scalingFactorsXrightAscensionContribution_ =
////            std::make_pair( - angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 0 ] * scalingFactorXCoefficient_.first
////                            + angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 1 ] * scalingFactorXCoefficient_.second,
////                            angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 0 ] * scalingFactorXCoefficient_.first
////                            + angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 1 ] * scalingFactorXCoefficient_.second );
//////    partialsXrightAscensionContribution =
//////            ( angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 ) - angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) )
//////            * std::cos( averageDeclination_ );

//////    scalingFactorsXdeclinationContribution_ =
//////            std::make_pair( - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
//////                            * std::sin( averageDeclination_ ) * angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ),
//////                            - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
//////                            * std::sin( averageDeclination_ ) * angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) );
////////    partialsXdeclinationContribution *=  - ( rightAscensionSecondTransmitter_- rightAscensionFirstTransmitter_ ) / 2.0
////////            * std::sin( averageDeclination_ );
////////            - ( rightAscensionAndDeclinationSecondTransmitter.first - rightAscensionAndDeclinationFirstTransmitter.first ) / 2.0
////////            * std::sin( averageDeclination_ )
////////            * ( angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) + angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) );

////    YscalingFactors_ =
////            std::make_pair( - angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorYCoefficient_,
////            angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) * scalingFactorYCoefficient_ );

////    YlightTimeCorrectionScalingFactors_ =
////            std::make_pair( - angularPositionLightTimeCorrectionScalingFirstTransmitter_[ 1 ] * scalingFactorYCoefficient_,
////            angularPositionLightTimeCorrectionScalingSecondTransmitter_[ 1 ] * scalingFactorYCoefficient_ );

////    //Eigen::Matrix< double, 1, 3 >::Zero( 1, 3 );
//////    partialsY = angularPositionScalingFactorSecondTransmitter_.block( 1, 0, 1, 3 ) - angularPositionScalingFactorFirstTransmitter_.block( 1, 0, 1, 3 );

//////    double X = std::cos( averageDeclination_ )
//////            * ( rightAscensionAndDeclinationSecondTransmitter.first - rightAscensionAndDeclinationFirstTransmitter.first );

//////    double Y = rightAscensionAndDeclinationSecondTransmitter.second - rightAscensionAndDeclinationFirstTransmitter.second;

//////    referenceScalingFactor_ = 1.0 / std::sqrt( X * X + Y * Y )
//////            * ( X * partialsXWrtLinkEndPosition + Y * partialsY );

//////    // UPDATE REFERENCE LIGHT TIME CORRECTION SCALING

    currentLinkEndType_ = fixedLinkEnd;

}

//! Function to calculate the observation partial(s) at required time and state
MutualApproximationPartial::MutualApproximationPartialReturnType MutualApproximationPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    if( linkEndOfFixedTime != mutualApproximationScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error mutual approximation partial and scaling are inconsistent" );
    }

    MutualApproximationPartialReturnType returnPartial;

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

//        // Scale position partials
//        returnPartial.push_back(
//                    std::make_pair(
//                        mutualApproximationScaler_->getScalingFactor( positionPartialIterator_->first ) *
//                        ( positionPartialIterator_->second->calculatePartialOfPosition(
//                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }


//    // Add scaled light-time correcion partials.
//    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ); i++ )
//    {
//        currentLinkTimeCorrectionPartialFirstTransmitter_ = lighTimeCorrectionPartialsFunctionsFirstTransmitter_.at( i )( states, times );
//        currentLinkTimeCorrectionPartialSecondTransmitter_ = lighTimeCorrectionPartialsFunctionsSecondTransmitter_.at( i )( states, times );

//        if ( currentLinkTimeCorrectionPartialFirstTransmitter_.second !=
//             currentLinkTimeCorrectionPartialSecondTransmitter_.second )
//        {
//            throw std::runtime_error( "Error when making mutual approximation light time correction partials, unconsistency"
//                                      " in receiver times between receiver - first transmitter and receiver - second transmitter legs." );
//        }

//        returnPartial.push_back(
//                    std::make_pair( mutualApproximationScaler_->getLightTimePartialScalingFactor( ) *
//                                    physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartialFirstTransmitter_.first,
//                    currentLinkTimeCorrectionPartialFirstTransmitter_.second ) );
//    }


    return returnPartial;
}

}

}
