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
        const Eigen::Vector3d& cartesianPositionVector )
{

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - cartesianPositionVector( 1 );
    partial( 1 ) = cartesianPositionVector( 0 );
    partial /=  ( cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) +
                                       cartesianPositionVector( 1 ) * cartesianPositionVector( 1 ) );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d cartesianPositionVector )
{
    // Set partial vector
    double range = cartesianPositionVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - cartesianPositionVector( 0 ) * cartesianPositionVector( 2 );
    partial( 1 ) = - cartesianPositionVector( 1 ) * cartesianPositionVector( 2 );
    partial( 2 ) = cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) + cartesianPositionVector( 1 ) * cartesianPositionVector( 1 );
    partial /= ( range * range * std::sqrt( cartesianPositionVector( 0 ) * cartesianPositionVector( 0 ) + cartesianPositionVector( 1 ) * cartesianPositionVector( 1 ) ) );

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


//! Function to compute the partial w.r.t. velocity of observer or observed object of the first time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndVelocity(
        const Eigen::Vector3d& cartesianPositionVector )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );

    Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
    partials( 0, 0 ) = - ry / ( rx * rx + ry * ry );
    partials( 0, 1 ) = rx / ( rx * rx + ry * ry );

    return partials;
}


//! Function to compute the partial w.r.t. velocity of observer or observed object of the first time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndVelocity(
        const Eigen::Vector3d& cartesianPositionVector )
{
    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );
    double rz = cartesianPositionVector.z( );

    Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
    partials( 0, 0 ) = - rx * rz;
    partials( 0, 1 ) = - ry * rz;
    partials( 0, 2 ) = rx * rx + ry * ry;

    return partials / ( std::sqrt( rx * rx + ry * ry ) * ( rx * rx + ry * ry + rz * rz ) );
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

    double partial = 2.0 * ( ( 1.0 / ( ( intermediateVariable * intermediateVariable + ry * ry ) * ( intermediateVariable * intermediateVariable + ry * ry ) ) )
                      * ( 1.0 / std::sqrt( rx * rx + ry * ry ) ) * ( rx * vy - ry * vx )
                      * ( ( ry * ry - intermediateVariable * intermediateVariable ) * vx - 2.0 * intermediateVariable * ry * vy
                          + ( 1.0 / ( rx * rx + ry * ry ) ) * ( rx * vx + ry * vy )
                          * ( intermediateVariable * intermediateVariable * ( - 2.0 * std::sqrt( rx * rx + ry * ry ) - rx ) - rx * ry * ry ) )
                      + intermediateVariable / ( intermediateVariable * intermediateVariable + ry * ry ) * 1.0 / std::sqrt( rx * rx + ry * ry )
                      * ( rx * ay - ry * ax ) );

    partial = - 2.0 * ( rx * vy - ry * vx ) * ( rx * vx + ry * vy ) / ( ( rx * rx + ry * ry ) * ( rx * rx + ry * ry ) )
            + ( rx * ay - ry * ax ) / ( rx * rx + ry * ry );

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

    partial = 1.0 / ( std::sqrt( rx * rx + ry * ry ) * ( rx * rx + ry * ry + rz * rz ) )
            * ( ( - rx * rz * ax - ry * rz * ay + ( rx * rx + ry * ry ) * az )
                - rz * ( rx * vy - ry * vx ) * ( rx * vy - ry * vx ) / ( rx * rx + ry * ry )
            + 2.0 * ( rx * vx + ry * vy + rz * vz ) / ( rx * rx + ry * ry + rz * rz )
                * ( rz * ( rx * vx + ry * vy ) - vz * ( rx * rx + ry * ry ) ) );

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
        partial( 0 ) = 2.0 * ( intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
                               * ( ay + rx * accelerationPartialWrtX.y( ) - ry * accelerationPartialWrtX.x( ) ) );

        partial( 1 ) = 2.0 * ( intermediateVariableA / ( intermediateVariableA * intermediateVariableA + ry * ry ) * 1.0 / intermediateVariableC
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
        partial( 0 ) = 1.0 / intermediateVariableC * accelerationPartialWrtX.z( )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ax + rx * accelerationPartialWrtX.x( ) + ry * accelerationPartialWrtX.y( )
                                                     + rz * accelerationPartialWrtX.z( ) );

        partial( 1 ) = 1.0 / intermediateVariableC * accelerationPartialWrtY.z( )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( ay + rx * accelerationPartialWrtY.x( ) + ry * accelerationPartialWrtY.y( )
                                                     + rz * accelerationPartialWrtY.z( ) );

        partial( 2 ) = 1.0 / intermediateVariableC * accelerationPartialWrtZ.z( )
                - rz / ( intermediateVariableC * intermediateVariableD ) * ( az + rx * accelerationPartialWrtZ.x( ) + ry * accelerationPartialWrtZ.y( )
                                                     + rz * accelerationPartialWrtZ.z( ) );
    }

    if ( !wrtTransmitterPosition )
    {
        partial *= - 1.0;
    }

    return partial;
}


//! Function to compute the partial w.r.t. velocity of observer or observed object of the second time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
        const Eigen::Matrix3d& partialCartesianAccelerationWrtVelocity,
        const bool wrtTransmitterPosition,
        const bool wrtOtherLinkEnd )
{

    double rx = cartesianPositionVector.x( );
    double ry = cartesianPositionVector.y( );

    double vx = cartesianVelocityVector.x( );
    double vy = cartesianVelocityVector.y( );

    Eigen::Vector3d accelerationPartialWrtVx = partialCartesianAccelerationWrtVelocity.block( 0, 0, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtVy = partialCartesianAccelerationWrtVelocity.block( 0, 1, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtVz = partialCartesianAccelerationWrtVelocity.block( 0, 2, 3, 1 );
    if ( !wrtTransmitterPosition )
    {
        accelerationPartialWrtVx *= - 1.0;
        accelerationPartialWrtVy *= - 1.0;
        accelerationPartialWrtVz *= - 1.0;
    }

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );

    if ( !wrtOtherLinkEnd )
    {
        partial( 0 ) = 2.0 / ( ( rx * rx + ry * ry ) * ( rx * rx + ry * ry ) )
                * ( 2.0 * rx * ry * vx + ( ry * ry - rx * rx ) * vy )
                + 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVx[ 1 ] - ry * accelerationPartialWrtVx[ 0 ] );


        partial( 1 ) = 2.0 / ( ( rx * rx + ry * ry ) * ( rx * rx + ry * ry ) )
                * ( - 2.0 * rx * ry * vy + ( ry * ry - rx * rx ) * vx )
                + 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVy[ 1 ] - ry * accelerationPartialWrtVy[ 0 ] );


        partial( 2 ) = 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVz[ 1 ] - ry * accelerationPartialWrtVz[ 0 ] );
    }

    else
    {
        partial( 0 ) = 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVx[ 1 ] - ry * accelerationPartialWrtVx[ 0 ] );

        partial( 1 ) = 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVy[ 1 ] - ry * accelerationPartialWrtVy[ 0 ] );

        partial( 2 ) = 1.0 / ( rx * rx + ry * ry ) * ( rx * accelerationPartialWrtVz[ 1 ] - ry * accelerationPartialWrtVz[ 0 ] );
    }


    if ( !wrtTransmitterPosition )
    {
        partial *= - 1.0;
    }

    return partial;

}


//! Function to compute the partial w.r.t. velocity of observer or observed object of the second time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
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

    Eigen::Vector3d accelerationPartialWrtVx = partialCartesianAccelerationWrtPosition.block( 0, 0, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtVy = partialCartesianAccelerationWrtPosition.block( 0, 1, 3, 1 );
    Eigen::Vector3d accelerationPartialWrtVz = partialCartesianAccelerationWrtPosition.block( 0, 2, 3, 1 );
    if ( !wrtTransmitterPosition )
    {
        accelerationPartialWrtVx *= - 1.0;
        accelerationPartialWrtVy *= - 1.0;
        accelerationPartialWrtVz *= - 1.0;
    }

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );


    if ( !wrtOtherLinkEnd )
    {
        partial( 0 ) = ( - rx * rz * accelerationPartialWrtVx[ 0 ] - ry * rz * accelerationPartialWrtVx[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVx[ 2 ] )
                + 2.0 * rz * ( ry * vx - rx * vy ) / ( rx * rx + ry * ry ) * ( - ry )
                + 2.0 / ( rx * rx + ry * ry + rz * rz )
                * ( - rx * vz * ( rx * rx + ry * ry - rz * rz )
                    + 2.0 * rx * rz * ( rx * vx + ry * vy ) );


        partial( 1 ) = ( - rx * rz * accelerationPartialWrtVy[ 0 ] - ry * rz * accelerationPartialWrtVy[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVy[ 2 ] )
                + 2.0 * rz * ( ry * vx - rx * vy ) / ( rx * rx + ry * ry ) * ( rx )
                + 2.0 / ( rx * rx + ry * ry + rz * rz )
                * ( - ry * vz * ( rx * rx + ry * ry - rz * rz )
                    + 2.0 * ry * rz * ( rx * vx + ry * vy ) );

        partial( 2 ) = ( - rx * rz * accelerationPartialWrtVz[ 0 ] - ry * rz * accelerationPartialWrtVz[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVz[ 2 ] )
                + 2.0 / ( rx * rx + ry * ry + rz * rz )
                * ( - ( rx * rx + ry * ry - rz * rz ) * ( rx * vx + ry * vy )
                    - 2.0 * rz * vz * ( rx * rx + ry * ry ) );
    }

    else
    {
        partial( 0 ) = ( - rx * rz * accelerationPartialWrtVx[ 0 ] - ry * rz * accelerationPartialWrtVx[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVx[ 2 ] );

        partial( 1 ) = ( - rx * rz * accelerationPartialWrtVy[ 0 ] - ry * rz * accelerationPartialWrtVy[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVy[ 2 ] );

        partial( 2 ) = ( - rx * rz * accelerationPartialWrtVz[ 0 ] - ry * rz * accelerationPartialWrtVz[ 1 ]
                + ( rx * rx + ry * ry ) * accelerationPartialWrtVz[ 2 ] );
    }

    partial /= ( std::sqrt( rx * rx + ry * ry ) * ( rx * rx + ry * ry + rz * rz ) );

    if ( !wrtTransmitterPosition )
    {
        partial *= - 1.0;
    }

    return partial;
}



//! Compute partial of instrumental frame relative position w.r.t. right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativePositionWrtRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // X w.r.t. first transmitter right ascension.
    partials( 0, 0 ) = - std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // X w.r.t. second transmitter right ascension.
    partials( 0, 1 ) = std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Y w.r.t. first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Y w.r.t. second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}


//! Compute partial of instrumental frame relative position w.r.t. declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativePositionWrtDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // X w.r.t. first transmitter declination.
    partials( 0, 0 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // X w.r.t. second transmitter declination.
    partials( 0, 1 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Y w.r.t. first transmitter declination.
    partials( 1, 0 ) = - 1.0;

    // Y w.r.t. second transmitter declination.
    partials( 1, 1 ) = 1.0;

    return partials;
}


//! Compute partial of instrumental frame relative velocity w.r.t. right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeVelocityWrtRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Vx w.r.t. first transmitter right ascension.
    partials( 0, 0 ) = 1.0 / 2.0 * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Vx w.r.t. second transmitter right ascension.
    partials( 0, 1 ) = - 1.0 / 2.0 * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Vy w.r.t. first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Vy w.r.t. second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}

////! Compute partial of instrumental frame relative velocity w.r.t. right ascension of second transmitter.
//Eigen::Matrix< double, 2, 1 > computePartialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterRightAscension( )
//{
//}

//! Compute partial of instrumental frame relative velocity w.r.t. declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeVelocityWrtDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Vx w.r.t. first transmitter declination.
    partials( 0, 0 ) = - 1.0 / 2.0 * ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Vx w.r.t. second transmitter declination.
    partials( 0, 1 ) = - 1.0 / 2.0 * ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Vy w.r.t. first transmitter declination.
    partials( 1, 0 ) = 0.0;

    // Vy w.r.t. second transmitter declination.
    partials( 1, 1 ) = 0.0;

    return partials;
}

////! Compute partial of instrumental frame relative velocity w.r.t. declination of second transmitter.
//Eigen::Matrix< double, 2, 1 > computePartialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterDeclination( )
//{

//}

//! Compute partial of instrumental frame relative velocity w.r.t. first time derivative of right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Vx w.r.t. first time derivative of first transmitter right ascension.
    partials( 0, 0 ) = - std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Vx w.r.t. first time derivative of second transmitter right ascension.
    partials( 0, 1 ) = std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Vy w.r.t. first time derivative of first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Vy w.r.t. first time derivative of second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}

////! Compute partial of instrumental frame relative velocity w.r.t. first time derivative of right ascension of second transmitter.
//Eigen::Matrix< double, 2, 1 > computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeSecondTransmitterRightAscension( )
//{

//}

//! Compute partial of instrumental frame relative velocity w.r.t. first time derivative of declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Vx w.r.t. first time derivative of first transmitter declination.
    partials( 0, 0 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Vx w.r.t. first time derivative of second transmitter declination.
    partials( 0, 1 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Vy w.r.t. first time derivative of first transmitter declination.
    partials( 1, 0 ) = - 1.0;

    // Vy w.r.t. first time derivative of second transmitter declination.
    partials( 1, 1 ) = 1.0;

    return partials;
}

////! Compute partial of instrumental frame relative velocity w.r.t. first time derivative of declination of second transmitter.
//Eigen::Matrix< double, 2, 1 > computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeSecondTransmitterDeclination( )
//{

//}


//! Compute partial of instrumental frame relative acceleration w.r.t. right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. first transmitter right ascension.
    partials( 0, 0 ) = ( 1.0 / 4.0 ) * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            + ( 1.0 / 2.0 ) * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( secondPartialOfDeclinationFirstTransmitterWrtTime_ + secondPartialOfDeclinationSecondTransmitterWrtTime_ );

    // Ax w.r.t. second transmitter right ascension.
    partials( 0, 1 ) = - ( 1.0 / 4.0 ) * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            - ( 1.0 / 2.0 ) * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( secondPartialOfDeclinationFirstTransmitterWrtTime_ + secondPartialOfDeclinationSecondTransmitterWrtTime_ );

    // Ay w.r.t. first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Ay w.r.t. second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}


//! Compute partial of instrumental frame relative acceleration w.r.t. declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. first transmitter declination.
    partials( 0, 0 ) = - ( 1.0 / 2.0 ) * (
                ( secondPartialOfRightAscensionSecondTransmitterWrtTime_ - secondPartialOfRightAscensionFirstTransmitterWrtTime_ )
                - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ ) )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            + ( 1.0 / 2.0 )
            * ( - ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
                - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
                * ( secondPartialOfDeclinationFirstTransmitterWrtTime_ + secondPartialOfDeclinationSecondTransmitterWrtTime_ ) )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ax w.r.t. second transmitter declination.
    partials( 0, 1 ) = - ( 1.0 / 2.0 ) * (
                ( secondPartialOfRightAscensionSecondTransmitterWrtTime_ - secondPartialOfRightAscensionFirstTransmitterWrtTime_ )
                - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 4.0
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ ) )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            + ( 1.0 / 2.0 )
            * ( - ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
                * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
                - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
                * ( secondPartialOfDeclinationFirstTransmitterWrtTime_ + secondPartialOfDeclinationSecondTransmitterWrtTime_ ) )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ay w.r.t. first transmitter declination.
    partials( 1, 0 ) = 0.0;

    // Ay w.r.t. second transmitter declination.
    partials( 1, 1 ) = 0.0;

    return partials;
}


//! Compute partial of instrumental frame relative acceleration w.r.t. first time derivative of right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. first time derivative of first transmitter right ascension.
    partials( 0, 0 ) = std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Ax w.r.t. first time derivative of second transmitter right ascension.
    partials( 0, 1 ) = - std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    // Ay w.r.t. first time derivative of first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Ay w.r.t. first time derivative of second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}

//! Compute partial of instrumental frame relative acceleration w.r.t. first time derivative of declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. first time derivative of first transmitter declination.
    partials( 0, 0 ) = - ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ax w.r.t. first time derivative of second transmitter declination.
    partials( 0, 1 ) = - ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
            * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ )
            * std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ay w.r.t. first time derivative of first transmitter declination.
    partials( 1, 0 ) = 0.0;

    // Ay w.r.t. first time derivative of second transmitter declination.
    partials( 1, 1 ) = 0.0;

    return partials;
}


//! Compute partial of instrumental frame relative acceleration w.r.t. second time derivative of right ascension of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeRightAscension( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. second time derivative of first transmitter right ascension.
    partials( 0, 0 ) = - std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ax w.r.t. second time derivative of second transmitter right ascension.
    partials( 0, 1 ) = std::cos( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ay w.r.t. second time derivative of first transmitter right ascension.
    partials( 1, 0 ) = 0.0;

    // Ay w.r.t. second time derivative of second transmitter right ascension.
    partials( 1, 1 ) = 0.0;

    return partials;
}

//! Compute partial of instrumental frame relative acceleration w.r.t. second time derivative of declination of the transmitters.
Eigen::Matrix< double, 2, 2 > MutualApproximationScalingBase::computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeDeclination( )
{
    Eigen::Matrix< double, 2, 2 > partials;

    // Ax w.r.t. second time derivative of first transmitter declination.
    partials( 0, 0 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ax w.r.t. second time derivative of second transmitter declination.
    partials( 0, 1 ) = - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( ( declinationFirstTransmitter_ + declinationSecondTransmitter_ ) / 2.0 );

    // Ay w.r.t. second time derivative of first transmitter declination.
    partials( 1, 0 ) = - 1.0;

    // Ay w.r.t. second time derivative of second transmitter declination.
    partials( 1, 1 ) = 1.0;

    return partials;
}


//! Compute the angle theta used to derive the real solutions of the depressed cubic equation.
double computeAngleThetaRealSolutionsCubicEquation( double intermediateQ,
                                                    double intermediateR )
{
    return std::acos( intermediateR / std::sqrt( - intermediateQ * intermediateQ * intermediateQ ) );
}



//! Compute the partial of the intermediate variable Q (used to find the roots of the central instant depressed cubic polynomial)
//! w.r.t. link end state (either position of velocity). The partials are returned by reference.
void computePartialOfIntermediateVariableQWrtLinkEndState(
        Eigen::Vector3d depressedCubicPolynomialCoefficients,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterState,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterState,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableQwrtReceiverState )
{
    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 3, 3 > partialsOfDepressedPolynomialCoefficientsWrtLinkEndState;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterState;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterState;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtReceiverState;
        }


        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );

        // Retrieve partials of the second order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndState
                = partialsOfDepressedPolynomialCoefficientsWrtLinkEndState.block( 0, 0, 1, 3 );

        // Retrieve partials of the first order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndState
                = partialsOfDepressedPolynomialCoefficientsWrtLinkEndState.block( 1, 0, 1, 3 );


        // Compute partials of intermediate variable Q w.r.t. link end state.
        partials = ( 1.0 / 9.0 ) * ( 3.0 * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndState
                                           - 2.0 * depressedCubicPolynomialCoefficients[ 0 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndState );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialOfIntermediateVariableQwrtFirstTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialOfIntermediateVariableQwrtSecondTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver state.
        {
            partialOfIntermediateVariableQwrtReceiverState = partials;
        }

    }

}



//! Compute the partial of the intermediate variable R (used to find the roots of the central instant depressed cubic polynomial)
//! w.r.t. link end state (either position or velocity).
void computePartialOfIntermediateVariableRWrtLinkEndState(
        Eigen::Vector3d depressedCubicPolynomialCoefficients,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterState,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterState,
        Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfIntermediateVariableRwrtReceiverState )
{

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 3, 3 > partialsOfDepressedPolynomialCoefficientsWrtLinkEndState;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterState;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterState;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver state.
        {
            partialsOfDepressedPolynomialCoefficientsWrtLinkEndState = partialOfDepressedPolynomialCoefficientsWrtReceiverState;
        }


        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );

        // Retrieve partials of the second order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndState =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndState.block( 0, 0, 1, 3 );

        // Retrieve partials of the first order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndState =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndState.block( 1, 0, 1, 3 );

        // Retrieve partials of the zero order coefficient of the depressed cubic polynomial.
        Eigen::Matrix< double, 1, 3 > partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndState =
                partialsOfDepressedPolynomialCoefficientsWrtLinkEndState.block( 2, 0, 1, 3 );


        // Compute partials of intermediate variable Q w.r.t. link end state.
        partials = ( 1.0 / 54.0 ) * ( 9.0 * ( depressedCubicPolynomialCoefficients[ 0 ] * partialsOfFirstOrderCoefficientDepressedPolynomialWrtLinkEndState
                                            + depressedCubicPolynomialCoefficients[ 1 ] * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndState )
                                           - 27.0 * partialsOfZeroOrderCoefficientDepressedPolynomialWrtLinkEndState
                - 6.0 * depressedCubicPolynomialCoefficients[ 0 ] * depressedCubicPolynomialCoefficients[ 0 ]
                * partialsOfSecondOrderCoefficientDepressedPolynomialWrtLinkEndState );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialOfIntermediateVariableRwrtFirstTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialOfIntermediateVariableRwrtSecondTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver state.
        {
            partialOfIntermediateVariableRwrtReceiverState = partials;
        }

    }

}


//! Compute the partial of the angle theta (used to derive the real solutions of the central instant depressed cubic equation)
//! w.r.t. link end state (either position or velocity).
void computePartialOfAngleThetaCubicEquationWrtLinkEndState(
        double intermediateVariableQ,
        double intermediateVariableR,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverState,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialOfAngleThetaCubicEquationWrtReceiverState )
{

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndState;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndState;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtFirstTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtFirstTransmitterState;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtSecondTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtSecondTransmitterState;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver state
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtReceiverState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtReceiverState;
        }


        // Compute partials of intermediate variable Q w.r.t. link end state, for the two link legs
        // (first transmitter - receiver and second transmitter - receiver)
        Eigen::Matrix< double, 1, 3 > partials = std::sqrt( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
                                    / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ + intermediateVariableR * intermediateVariableR ) )
                * ( 1.0 / ( intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) )
                * ( std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) * partialsOfIntermediateVariableRwrtLinkEndState
                    + 3.0 * intermediateVariableR * intermediateVariableQ * intermediateVariableQ
                    / ( 2.0 * std::sqrt( - intermediateVariableQ * intermediateVariableQ * intermediateVariableQ ) ) * partialsOfIntermediateVariableQwrtLinkEndState );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialOfAngleThetaCubicEquationWrtFirstTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialOfAngleThetaCubicEquationWrtSecondTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver state.
        {
            partialOfAngleThetaCubicEquationWrtReceiverState = partials;
        }

    }

}


//! Compute the partial of the intermediate variable T (used to find the roots of the central instant depressed cubic polynomial)
//! w.r.t. link end state (either position or velocity).
void computePartialOfIntermediateVariableTWrtLinkEndState(
        double intermediateVariableQ,
        double intermediateVariableR,
        double intermediateVariableT,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableTwrtReceiverState )
{
    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndState;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndState;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtFirstTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtFirstTransmitterState;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtSecondTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtSecondTransmitterState;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtReceiverState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtReceiverState;
        }

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableBetaWrtLinkEndState =
                3.0 * intermediateVariableQ * intermediateVariableQ * partialsOfIntermediateVariableQwrtLinkEndState
                + 2.0 * intermediateVariableR * partialsOfIntermediateVariableRwrtLinkEndState;


        // Compute partials of intermediate variable Q w.r.t. link end state.
        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
        if ( intermediateVariableT >= 0.0 )
        {
            partials = 1.0 / ( 3.0 * std::pow( intermediateVariableR - std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndState - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndState );
        }
        else
        {
            partials = 1.0 / ( 3.0 * std::pow( - ( intermediateVariableR - std::sqrt( intermediateVariableBeta ) ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndState - 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndState );
        }


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfIntermediateVariableTwrtFirstTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfIntermediateVariableTwrtSecondTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver state.
        {
            partialsOfIntermediateVariableTwrtReceiverState = partials;
        }
    }

}



//! Compute the partial of the intermediate variable S (used to find the roots of the central instant depressed cubic polynomial)
//! w.r.t. link end state (either position or velocity).
void computePartialOfIntermediateVariableSWrtLinkEndState(
        double intermediateVariableQ,
        double intermediateVariableR,
        double intermediateVariableS,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtReceiverState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtReceiverState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtFirstTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtSecondTransmitterState,
        Eigen::Matrix< double, 1, 3 >& partialsOfIntermediateVariableSwrtReceiverState )
{
    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableQwrtLinkEndState;
        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableRwrtLinkEndState;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtFirstTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtFirstTransmitterState;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtSecondTransmitterState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtSecondTransmitterState;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver state.
        {
            partialsOfIntermediateVariableQwrtLinkEndState = partialsOfIntermediateVariableQwrtReceiverState;
            partialsOfIntermediateVariableRwrtLinkEndState = partialsOfIntermediateVariableRwrtReceiverState;
        }

        Eigen::Matrix< double, 1, 3 > partialsOfIntermediateVariableBetaWrtLinkEndState = 3.0 * intermediateVariableQ * intermediateVariableQ
                * partialsOfIntermediateVariableQwrtLinkEndState + 2.0 * intermediateVariableR * partialsOfIntermediateVariableRwrtLinkEndState;

        // Compute partials of intermediate variable Q w.r.t. link end state.
        Eigen::Matrix< double, 1, 3 > partials = Eigen::Matrix< double, 1, 3 >::Zero( );
        if ( intermediateVariableS >= 0.0 )
        {
            partials = 1.0 / ( 3.0 * std::pow( intermediateVariableR + std::sqrt( intermediateVariableBeta ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndState + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndState );
        }
        else
        {
            partials = 1.0 / ( 3.0 * std::pow( - ( intermediateVariableR + std::sqrt( intermediateVariableBeta ) ), 2.0 / 3.0 ) )
                    * ( partialsOfIntermediateVariableRwrtLinkEndState + 1.0 / ( 2.0 * std::sqrt( intermediateVariableBeta ) )
                        * partialsOfIntermediateVariableBetaWrtLinkEndState );
        }


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter state.
        {
            partialsOfIntermediateVariableSwrtFirstTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter state.
        {
            partialsOfIntermediateVariableSwrtSecondTransmitterState = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver state.
        {
            partialsOfIntermediateVariableSwrtReceiverState = partials;
        }
    }

}


Eigen::Vector2d MutualApproximationScalingBase::computeRelativePositionInInstrumentalFrame( )
{
    Eigen::Vector2d relativePositionInInstrumentalFrame = Eigen::Vector2d::Zero( );

    relativePositionInInstrumentalFrame[ 0 ] = ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
            * std::cos( averageDeclination_ );
    relativePositionInInstrumentalFrame[ 1 ] = declinationSecondTransmitter_ - declinationFirstTransmitter_;

    return relativePositionInInstrumentalFrame;
}


Eigen::Vector2d MutualApproximationScalingBase::computeRelativeVelocityInInstrumentalFrame( )
{
    Eigen::Vector2d relativeVelocityInInstrumentalFrame = Eigen::Vector2d::Zero( );

    relativeVelocityInInstrumentalFrame[ 0 ] = ( partialOfRightAscensionSecondTransmitterWrtTime_ - partialOfRightAscensionFirstTransmitterWrtTime_ )
            * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0
            * std::sin( averageDeclination_ ) * ( partialOfDeclinationFirstTransmitterWrtTime_ + partialOfDeclinationSecondTransmitterWrtTime_ );

    relativeVelocityInInstrumentalFrame[ 1 ] = partialOfDeclinationSecondTransmitterWrtTime_ - partialOfDeclinationFirstTransmitterWrtTime_;

    return relativeVelocityInInstrumentalFrame;
}


Eigen::Vector2d MutualApproximationScalingBase::computeRelativeAccelerationInInstrumentalFrame( )
{

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



void MutualApproximationScalingBase::computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPosition( )
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
    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 ) =
            - partialOfRightAscensionWrtPositionFirstTransmitter * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfDeclinationWrtPositionFirstTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 ) = - partialOfDeclinationWrtPositionFirstTransmitter;

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
    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 ) =
            partialOfRightAscensionWrtPositionSecondTransmitter * std::cos( averageDeclination_ )
            - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ ) / 2.0 * std::sin( averageDeclination_ )
            * partialOfDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 ) = partialOfDeclinationWrtPositionSecondTransmitter;

    // Compute contribution of the link second transmitter-receiver to the partials of relative position
    // (X and Y coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_ -= partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;


}


void MutualApproximationScalingBase::computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPosition( )
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
    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 ) =
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

    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 ) = - partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

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
    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 ) =
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

    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 ) = partialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition;

    // Compute contribution of the link second transmitter-receiver to the partials of relative velocity
    // (Vx and Vy coordinates in the instrumental frame of the receiver) w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_ -= partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;

}


void MutualApproximationScalingBase::computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
        Eigen::Vector3d relativeAccelerationFirstTransmitterWrtReceiver,
        Eigen::Vector3d relativeAccelerationSecondTransmitterWrtReceiver,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterPosition )
{

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
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 ) =
            ( partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition - partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition ) * std::cos( averageDeclination_ )
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
            * ( partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition + partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition) ;

    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 ) = partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition
            - partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition;



    /// SAME W.R.T. RECEIVER POSITION FOR FIRST TRANSMITTER - RECEIVER LINK
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 0, 0, 1, 3 ) = - 1.0 * (
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

    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 1, 0, 1, 3 ) = - 1.0 * ( - partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition );




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
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition(
                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,                                                     relativeAccelerationFirstTransmitterWrtReceiver, partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );


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
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition(
                relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver, relativeAccelerationFirstTransmitterWrtReceiver,
                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );



    // Compute partials of the relative acceleration (Ax and Ay coordinates in the instrumental frame of the receiver)
    // w.r.t. link end position, for the receiver - second transmitter link.

    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 ) =
            ( partialOfSecondTimeDerivativeRightAscensionWrtTransmitterPosition - partialOfSecondTimeDerivativeRightAscensionWrtOtherTransmitterPosition ) * std::cos( averageDeclination_ )
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
            * ( partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition + partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition );

    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 ) =
            partialOfSecondTimeDerivativeDeclinationWrtTransmitterPosition - partialOfSecondTimeDerivativeDeclinationWrtOtherTransmitterPosition;


    /// SAME W.R.T. RECEIVER POSITION FOR SECOND TRANSMITTER - RECEIVER LINK (ADDED TO CONTRIBUTION OF THE FIRST TRANSMITTER - RECEIVER LINK)
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 0, 0, 1, 3 ) += - 1.0 * (
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

    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 1, 0, 1, 3 ) += - 1.0 * ( partialOfSecondTimeDerivativeDeclinationWrtReceiverPosition );

}


void MutualApproximationScalingBase::computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndVelocity( )
{

    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterVelocity_ = Eigen::Matrix< double, 2, 3 >::Zero( );
    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterVelocity_ = Eigen::Matrix< double, 2, 3 >::Zero( );
    partialsOfInstrumentalFrameRelativePositionWrtReceiverVelocity_ = Eigen::Matrix< double, 2, 3 >::Zero( );
}


void MutualApproximationScalingBase::computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndVelocity( )
{
    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

//    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
//    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeRightAscension( );

    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeDeclination( );

    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtVelocityFirstTransmitter
            = computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtVelocityFirstTransmitter
            = computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtVelocitySecondTransmitter
            = computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtVelocitySecondTransmitter
            = computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver );

    /// Partials w.r.t. first transmitter velocity.
    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity_.block( 0, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocityFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocityFirstTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity_.block( 1, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocityFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocityFirstTransmitter;

    /// Partials w.r.t. second transmitter velocity.
    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity_.block( 0, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocitySecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocitySecondTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity_.block( 1, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocitySecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocitySecondTransmitter;

    /// Partials w.r.t. receiver velocity.
    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity_.block( 0, 0, 1, 3 ) =
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocityFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocityFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocitySecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocitySecondTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity_.block( 1, 0, 1, 3 ) =
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocityFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocityFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtVelocitySecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtVelocitySecondTransmitter;

}


void MutualApproximationScalingBase::computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndVelocity(
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverVelocity,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterVelocity,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverVelocity,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterVelocity,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity )
{

    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeDeclination( );

    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeDeclination( );


    /// First transmitter - receiver link.

    // Compute partials of first time derivative of the first transmitter right ascension w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterVelocity =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver );

    // Compute partials of first time derivative of the first transmitter declination w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterVelocity =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver );


    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtTransmitterVelocity, true );

    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity, true, true );

    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. receiver velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtReceiverVelocity, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver velocity while the others are all w.r.t. transmitter velocity.
//    // This factor is needed when computing the velocity partial of the instrumental frame relative acceleration w.r.t. receiver velocity as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverVelocity *= - 1.0;


    // Compute partials of second time derivative of the first transmitter declination w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtTransmitterVelocity, true );

    // Compute partials of second time derivative of the second transmitter declination w.r.t. first transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity, true, true );

    // Compute partials of second time derivative of the first transmitter declination w.r.t. receiver velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtReceiverVelocity, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver velocity while the others are all w.r.t. transmitter velocity.
//    // This factor is needed when computing the position velocity of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverVelocity *= - 1.0;



    /// Second transmitter - receiver link.

    // Compute partials of first time derivative of the second transmitter right ascension w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterVelocity =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver );

    // Compute partials of first time derivative of the second transmitter declination w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterVelocity =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver );


    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtTransmitterVelocity, true );

    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity, true, true );

    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. receiver velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverVelocity =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtReceiverVelocity, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver velocity while the others are all w.r.t. transmitter velocity.
//    // This factor is needed when computing the velocity partial of the instrumental frame relative acceleration w.r.t. receiver velocity as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverVelocity *= - 1.0;


    // Compute partials of second time derivative of the second transmitter declination w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtTransmitterVelocity, true );

    // Compute partials of second time derivative of the first transmitter declination w.r.t. second transmitter velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionFirstTransmitterWrtReceiver, relativeVelocityFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity, true, true );

    // Compute partials of second time derivative of the second transmitter declination w.r.t. receiver velocity.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverVelocity =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndVelocity( relativePositionSecondTransmitterWrtReceiver, relativeVelocitySecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtReceiverVelocity, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver velocity while the others are all w.r.t. transmitter velocity.
//    // This factor is needed when computing the position velocity of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverVelocity *= - 1.0;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. first transmitter velocity
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterVelocity_.block( 0, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 0 )
            * partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 0 )
            * partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterVelocity;


    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterVelocity_.block( 1, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 0 )
            * partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 0 )
            * partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterVelocity;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. second transmitter velocity
    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterVelocity_.block( 0, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 1 )
            * partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 1 )
            * partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterVelocity;


    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterVelocity_.block( 1, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 1 )
            * partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 1 )
            * partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterVelocity;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. receiver velocity
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverVelocity_.block( 0, 0, 1, 3 )
            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 0 )
            * ( - partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 0 )
            * ( - partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 1 )
            * ( - partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 1 )
            * ( - partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterVelocity )


            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverVelocity;


    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverVelocity_.block( 1, 0, 1, 3 )
            = partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 0 )
            * ( - partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 0 )
            * ( - partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 1 )
            * ( - partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterVelocity )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 1 )
            * ( - partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterVelocity )


            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverVelocity

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverVelocity;

}


void MutualApproximationScalingBase::computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPositionV2( )
{
    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    // Partials of relative position in the instrumental frame w.r.t. right ascensions/declinations.
    Eigen::Matrix< double, 2, 2 > partialsOfRelativePositionWrtRightAscension
            = computePartialsOfInstrumentalFrameRelativePositionWrtRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativePositionWrtDeclination
            = computePartialsOfInstrumentalFrameRelativePositionWrtDeclination( );

    // Partials of right ascensions and declinations w.r.t. link end positions.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );


    /// Partials w.r.t. first transmitter position.
    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 ) =
            partialsOfRelativePositionWrtRightAscension( 0, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativePositionWrtDeclination( 0, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 ) =
            partialsOfRelativePositionWrtRightAscension( 1, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativePositionWrtDeclination( 1, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter;

    /// Partials w.r.t. second transmitter position.
    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 ) =
            partialsOfRelativePositionWrtRightAscension( 0, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativePositionWrtDeclination( 0, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 ) =
            partialsOfRelativePositionWrtRightAscension( 1, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativePositionWrtDeclination( 1, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter;

    /// Partials w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_.block( 0, 0, 1, 3 ) =
            - partialsOfRelativePositionWrtRightAscension( 0, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativePositionWrtDeclination( 0, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            - partialsOfRelativePositionWrtRightAscension( 0, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativePositionWrtDeclination( 0, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_.block( 1, 0, 1, 3 ) =
            - partialsOfRelativePositionWrtRightAscension( 1, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativePositionWrtDeclination( 1, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            - partialsOfRelativePositionWrtRightAscension( 1, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativePositionWrtDeclination( 1, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter;

}


void MutualApproximationScalingBase::computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPositionV2( )
{
    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );

    // Partials of relative velocity in the instrumental frame w.r.t. right ascensions/declinations and their time derivatives.
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtRightAscension
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtDeclination
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtDeclination( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeVelocityWrtFirstTimeDerivativeDeclination( );

    // Partials of right ascensions and declinations w.r.t. link end positions.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionFirstTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionFirstTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionWrtPositionSecondTransmitter =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationWrtPositionSecondTransmitter =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    // Partials of first time derivatives of right ascensions and declinations w.r.t. link end positions.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtPositionFirstTransmitter
            = computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                   relativeVelocityFirstTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtPositionFirstTransmitter
            = computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                relativeVelocityFirstTransmitterWrtReceiver );

    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionWrtPositionSecondTransmitter
            = computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                   relativeVelocitySecondTransmitterWrtReceiver );
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationWrtPositionSecondTransmitter
            = computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                relativeVelocitySecondTransmitterWrtReceiver );


    /// Partials w.r.t. first transmitter position.
    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtRightAscension( 0, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtDeclination( 0, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionFirstTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtRightAscension( 1, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtDeclination( 1, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionFirstTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionFirstTransmitter;

    /// Partials w.r.t. second transmitter position.
    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtRightAscension( 0, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtDeclination( 0, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 ) =
            partialsOfRelativeVelocityWrtRightAscension( 1, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtDeclination( 1, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionSecondTransmitter
            + partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionSecondTransmitter;

    /// Partials w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_.block( 0, 0, 1, 3 ) =
            - partialsOfRelativeVelocityWrtRightAscension( 0, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtDeclination( 0, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionFirstTransmitter

            - partialsOfRelativeVelocityWrtRightAscension( 0, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtDeclination( 0, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 0, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 0, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionSecondTransmitter;

    partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_.block( 1, 0, 1, 3 ) =
            - partialsOfRelativeVelocityWrtRightAscension( 1, 0 ) * partialOfRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtDeclination( 1, 0 ) * partialOfDeclinationWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 0 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionFirstTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 0 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionFirstTransmitter

            - partialsOfRelativeVelocityWrtRightAscension( 1, 1 ) * partialOfRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtDeclination( 1, 1 ) * partialOfDeclinationWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeRightAscension( 1, 1 ) * partialOfFirstTimeDerivativeRightAscensionWrtPositionSecondTransmitter
            - partialsOfRelativeVelocityWrtFirstTimeDerivativeDeclination( 1, 1 ) * partialOfFirstTimeDerivativeDeclinationWrtPositionSecondTransmitter;

}



void MutualApproximationScalingBase::computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPositionV2(
        Eigen::Vector3d relativeAccelerationFirstTransmitterWrtReceiver,
        Eigen::Vector3d relativeAccelerationSecondTransmitterWrtReceiver,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterPosition,
        Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterPosition,
        Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterPosition )
{

    Eigen::Vector3d relativePositionFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 );
    Eigen::Vector3d relativePositionSecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 );

    Eigen::Vector3d relativeVelocityFirstTransmitterWrtReceiver = ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 );
    Eigen::Vector3d relativeVelocitySecondTransmitterWrtReceiver = ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 );


    // Partials of relative velocity in the instrumental frame w.r.t. right ascensions/declinations and their time derivatives.
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtRightAscension
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtDeclination
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtDeclination( );

    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTimeDerivativeDeclination( );

    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeRightAscension( );
    Eigen::Matrix< double, 2, 2 > partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination
            = computePartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTimeDerivativeDeclination( );


    /// First transmitter - receiver link.

    // Compute partials of the first transmitter right ascension w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionFirstTransmitterWrtTransmitterPosition =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );

    // Compute partials of the first transmitter declination w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationFirstTransmitterWrtTransmitterPosition =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver );


    // Compute partials of first time derivative of the first transmitter right ascension w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                 relativeVelocityFirstTransmitterWrtReceiver );

    // Compute partials of first time derivative of the first transmitter declination w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                              relativeVelocityFirstTransmitterWrtReceiver );


    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                  relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtTransmitterPosition, true );

    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                  relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtOtherTransmitterPosition, true, true );

    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                  relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtReceiverPosition, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver position while the others are all w.r.t. transmitter position.
//    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver position as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverPosition *= - 1.0;


    // Compute partials of second time derivative of the first transmitter declination w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                               relativeVelocityFirstTransmitterWrtReceiver,
                                                                               relativeAccelerationFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtTransmitterPosition, true );

    // Compute partials of second time derivative of the second transmitter declination w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                               relativeVelocitySecondTransmitterWrtReceiver,
                                                                               relativeAccelerationSecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtOtherTransmitterPosition, true, true );

    // Compute partials of second time derivative of the first transmitter declination w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                               relativeVelocityFirstTransmitterWrtReceiver,
                                                                               relativeAccelerationFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtReceiverPosition, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver position while the others are all w.r.t. transmitter position.
//    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter position.
//    partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverPosition *= - 1.0;



    /// Second transmitter - receiver link.

    // Compute partials of the first transmitter right ascension w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfRightAscensionSecondTransmitterWrtTransmitterPosition =
            computePartialOfRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );

    // Compute partials of the first transmitter declination w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfDeclinationSecondTransmitterWrtTransmitterPosition =
            computePartialOfDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver );


    // Compute partials of first time derivative of the second transmitter right ascension w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterPosition =
            computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                 relativeVelocitySecondTransmitterWrtReceiver );

    // Compute partials of first time derivative of the second transmitter declination w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterPosition =
            computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                              relativeVelocitySecondTransmitterWrtReceiver );


    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                  relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtTransmitterPosition, true );

    // Compute partials of second time derivative of the first transmitter right ascension w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                                  relativeVelocityFirstTransmitterWrtReceiver,
                                                                                  relativeAccelerationFirstTransmitterWrtReceiver,
                                                                                  partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );

    // Compute partials of second time derivative of the second transmitter right ascension w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                                  relativeVelocitySecondTransmitterWrtReceiver,
                                                                                  relativeAccelerationSecondTransmitterWrtReceiver,
                                                                                  partialAccelerationSecondTransmitterWrtReceiverPosition, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver position while the others are all w.r.t. transmitter position.
//    // This factor is needed when computing the position partial of the instrumental frame relative acceleration w.r.t. receiver position as minus the partial w.r.t. transmitter velocity.
//    partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverPosition *= - 1.0;


    // Compute partials of second time derivative of the second transmitter declination w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                               relativeVelocitySecondTransmitterWrtReceiver,
                                                                               relativeAccelerationSecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtTransmitterPosition, true );

    // Compute partials of second time derivative of the first transmitter declination w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionFirstTransmitterWrtReceiver,
                                                                               relativeVelocityFirstTransmitterWrtReceiver,
                                                                               relativeAccelerationFirstTransmitterWrtReceiver,
                                                                               partialAccelerationFirstTransmitterWrtOtherTransmitterPosition, true, true );

    // Compute partials of second time derivative of the second transmitter declination w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverPosition =
            computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition( relativePositionSecondTransmitterWrtReceiver,
                                                                               relativeVelocitySecondTransmitterWrtReceiver,
                                                                               relativeAccelerationSecondTransmitterWrtReceiver,
                                                                               partialAccelerationSecondTransmitterWrtReceiverPosition, false );

//    // Multiply by factor -1 because this partial specifically is w.r.t. receiver position while the others are all w.r.t. transmitter position.
//    // This factor is needed when computing the position position of the instrumental frame relative acceleration w.r.t. receiver state as minus the partial w.r.t. transmitter position.
//    partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverPosition *= - 1.0;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. first transmitter position.
    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_.block( 0, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtRightAscension( 0, 0 )
            * partialOfRightAscensionFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtDeclination( 0, 0 )
            * partialOfDeclinationFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 0 )
            * partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 0 )
            * partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterPosition;


    partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_.block( 1, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtRightAscension( 1, 0 )
            * partialOfRightAscensionFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtDeclination( 1, 0 )
            * partialOfDeclinationFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 0 )
            * partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 0 )
            * partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtFirstTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtFirstTransmitterPosition;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. second transmitter position.
    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_.block( 0, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtRightAscension( 0, 1 )
            * partialOfRightAscensionSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtDeclination( 0, 1 )
            * partialOfDeclinationSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 1 )
            * partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 1 )
            * partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterPosition;


    partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_.block( 1, 0, 1, 3 )

            = partialsOfRelativeAccelerationWrtRightAscension( 1, 1 )
            * partialOfRightAscensionSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtDeclination( 1, 1 )
            * partialOfDeclinationSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 1 )
            * partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 1 )
            * partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtSecondTransmitterPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtSecondTransmitterPosition;



    /// Compute partials of relative acceleration in the instrumental frame w.r.t. receiver position.
    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 0, 0, 1, 3 )
            = partialsOfRelativeAccelerationWrtRightAscension( 0, 0 )
            * ( - partialOfRightAscensionFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtDeclination( 0, 0 )
            * ( - partialOfDeclinationFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtRightAscension( 0, 1 )
            * ( - partialOfRightAscensionSecondTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtDeclination( 0, 1 )
            * ( - partialOfDeclinationSecondTransmitterWrtTransmitterPosition )


            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 0 )
            * ( - partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 0 )
            * ( - partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 0, 1 )
            * ( - partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 0, 1 )
            * ( - partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterPosition )


            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 0, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 0, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverPosition;


    partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_.block( 1, 0, 1, 3 ) =

            partialsOfRelativeAccelerationWrtRightAscension( 1, 0 )
            * ( - partialOfRightAscensionFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtDeclination( 1, 0 )
            * ( - partialOfDeclinationFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtRightAscension( 1, 1 )
            * ( - partialOfRightAscensionSecondTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtDeclination( 1, 1 )
            * ( - partialOfDeclinationSecondTransmitterWrtTransmitterPosition )


            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 0 )
            * ( - partialOfFirstTimeDerivativeRightAscensionFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 0 )
            * ( - partialOfFirstTimeDerivativeDeclinationFirstTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeRightAscension( 1, 1 )
            * ( - partialOfFirstTimeDerivativeRightAscensionSecondTransmitterWrtTransmitterPosition )

            + partialsOfRelativeAccelerationWrtFirstTimeDerivativeDeclination( 1, 1 )
            * ( - partialOfFirstTimeDerivativeDeclinationSecondTransmitterWrtTransmitterPosition )


            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 0 )
            * partialOfSecondTimeDerivativeRightAscensionFirstTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeRightAscension( 1, 1 )
            * partialOfSecondTimeDerivativeRightAscensionSecondTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 0 )
            * partialOfSecondTimeDerivativeDeclinationFirstTransmitterWrtReceiverPosition

            + partialsOfRelativeAccelerationWrtSecondTimeDerivativeDeclination( 1, 1 )
            * partialOfSecondTimeDerivativeDeclinationSecondTransmitterWrtReceiverPosition;

}



Eigen::Vector4d MutualApproximationScalingBase::computeCubicPolynomialCoefficients( )
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


void MutualApproximationScalingBase::computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition(
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

            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_;

            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;
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



void MutualApproximationScalingBase::computePartialOfCubicPolynomialCoefficientsWrtCartesianVelocity(
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity,
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity,
        Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtReceiverVelocity )
{

    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity;
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity;
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtLinkEndVelocity;

    // Successively compute partials w.r.t. velocity of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        if ( currentLinkEndIndex == 0 ) // Partials w.r.t. first transmitter velocity.
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterVelocity_;
            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity_;
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 1 ) // Partials w.r.t. second transmitter velocity.
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterVelocity_;
            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity_;
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 2 ) // Partials w.r.t. receiver velocity.
        {
            partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverVelocity_;
            partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity_;
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtReceiverVelocity_;
        }


        // Compute partials.

        Eigen::Vector3d partialsThirdOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                  instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 0 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 0 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 1 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 1 ),

                instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 2 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsSecondOrderCoefficient =
                3.0 * ( Eigen::Vector3d( ) <<
                        instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 0 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 0 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 0 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 0 ),

                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 1 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 1 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 1 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 1 ),

                instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 2 )
                + instrumentalFrameRelativeAcceleration_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 2 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 2 )
                + instrumentalFrameRelativeAcceleration_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsFirstOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 0 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 0 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 0 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 0 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 1 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 1 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 1 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 1 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 0, 2 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeAccelerationWrtLinkEndVelocity( 1, 2 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 2 )
                + 2.0 * instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 2 ) ).finished( ); //.transpose( );


        Eigen::Vector3d partialsZeroOrderCoefficient =
                2.0 * ( Eigen::Vector3d( ) <<
                        instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 0 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 0 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 1 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 1 ),

                instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 0, 2 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativeVelocityWrtLinkEndVelocity( 1, 2 ) ).finished( ); //.transpose( );



        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtReceiverVelocity.block( 0, 0, 1, 3 ) = partialsThirdOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverVelocity.block( 1, 0, 1, 3 ) = partialsSecondOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverVelocity.block( 2, 0, 1, 3 ) = partialsFirstOrderCoefficient.transpose( );
            partialOfCubicPolynomialCoefficientsWrtReceiverVelocity.block( 3, 0, 1, 3 ) = partialsZeroOrderCoefficient.transpose( );
        }
    }
}



Eigen::Vector3d MutualApproximationScalingBase::computeDepressedCubicPolynomialCoefficients(  )
{
    return ( Eigen::Vector3d( ) <<
             cubicPolynomialCoefficients_[ 1 ] / cubicPolynomialCoefficients_[ 0 ],
            cubicPolynomialCoefficients_[ 2 ] / cubicPolynomialCoefficients_[ 0 ],
            cubicPolynomialCoefficients_[ 3 ] / cubicPolynomialCoefficients_[ 0 ] ).finished( );
}


void MutualApproximationScalingBase::computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianPosition(
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtReceiverPosition )
{
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtReceiverPosition;
    computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition( partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
                                                                     partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
                                                                     partialOfCubicPolynomialCoefficientsWrtReceiverPosition );

//    std::cout << "partials of cubic polynomial coefficients w.r.t. first transmitter position: "
//              << partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition.transpose( ) << "\n\n";
//    std::cout << "partials of cubic polynomial coefficients w.r.t. second transmitter position: "
//              << partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition.transpose( ) << "\n\n";
//    std::cout << "partials of cubic polynomial coefficients w.r.t. receiver transmitter position: "
//              << partialOfCubicPolynomialCoefficientsWrtReceiverPosition.transpose( ) << "\n\n";

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

        // Compute partials of the first order coefficient of the depressed cubic polynomial.
        partials.block( 1, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 2, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 2 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 ) );

        // Compute partials of the zero order coefficient of the depressed cubic polynomial.
        partials.block( 2, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 3, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 3 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndPosition.block( 0, 0, 1, 3 ) );


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
}


void MutualApproximationScalingBase::computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianVelocity(
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterVelocity,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterVelocity,
        Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtReceiverVelocity )
{
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity;
    Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtReceiverVelocity;
    computePartialOfCubicPolynomialCoefficientsWrtCartesianVelocity( partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity,
                                                                     partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity,
                                                                     partialOfCubicPolynomialCoefficientsWrtReceiverVelocity );

//    std::cout << "partials of cubic polynomial coefficients w.r.t. first transmitter velocity: "
//              << partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity.transpose( ) << "\n\n";
//    std::cout << "partials of cubic polynomial coefficients w.r.t. second transmitter velocity: "
//              << partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity.transpose( ) << "\n\n";
//    std::cout << "partials of cubic polynomial coefficients w.r.t. receiver transmitter velocity: "
//              << partialOfCubicPolynomialCoefficientsWrtReceiverVelocity.transpose( ) << "\n\n";

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 4, 3 > partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity = partialOfCubicPolynomialCoefficientsWrtFirstTransmitterVelocity;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity = partialOfCubicPolynomialCoefficientsWrtSecondTransmitterVelocity;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver velocity.
        {
            partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity = partialOfCubicPolynomialCoefficientsWrtReceiverVelocity;
        }

        Eigen::Matrix< double, 3, 3 > partials = Eigen::Matrix< double, 3, 3 >::Zero( );

        // Compute partials of the second order coefficient of the depressed cubic polynomial.
        partials.block( 0, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 1, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 1 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 0, 0, 1, 3 ) );

        // Compute partials of the first order coefficient of the depressed cubic polynomial.
        partials.block( 1, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 2, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 2 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 0, 0, 1, 3 ) );

        // Compute partials of the zero order coefficient of the depressed cubic polynomial.
        partials.block( 2, 0, 1, 3 ) =
                ( 1.0 / ( cubicPolynomialCoefficients_[ 0 ] * cubicPolynomialCoefficients_[ 0 ] ) )
                * ( cubicPolynomialCoefficients_[ 0 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 3, 0, 1, 3 )
                - cubicPolynomialCoefficients_[ 3 ] * partialOfCubicPolynomialCoefficientsWrtLinkEndVelocity.block( 0, 0, 1, 3 ) );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity.
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterVelocity = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity.
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterVelocity = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver velocity.
        {
            partialOfDepressedCubicPolynomialCoefficientsWrtReceiverVelocity = partials;
        }

    }
}


void MutualApproximationScalingBase::computePartialsOfCentralInstantWrtLinkEndPosition(
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
    computePartialOfIntermediateVariableQWrtLinkEndState(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
                partialsIntermediateVariableQwrtFirstTransmitterPosition, partialsIntermediateVariableQwrtSecondTransmitterPosition,
                partialsIntermediateVariableQwrtReceiverPosition );

    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtFirstTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtSecondTransmitterPosition;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtReceiverPosition;
    computePartialOfIntermediateVariableRWrtLinkEndState(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterPosition,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterPosition, partialOfDepressedPolynomialCoefficientsWrtReceiverPosition,
                partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
                partialsIntermediateVariableRwrtReceiverPosition );

    if ( intermediateVariableBeta < 0 )
    {
//        std::cout << "BETA NEGATITVE" << "\n\n";

        double thetaAngle = computeAngleThetaRealSolutionsCubicEquation( intermediateVariableQ, intermediateVariableR );

        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtFirstTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtSecondTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtReceiverPosition;
        computePartialOfAngleThetaCubicEquationWrtLinkEndState( intermediateVariableQ, intermediateVariableR,
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

//        std::cout << "estimated central instant: " << times[ 0 ] << "\n\n";
//        std::cout << "firstSolutionCentralInstant: " << firstSolutionCentralInstant << "\n\n";
//        std::cout << "secondSolutionCentralInstant: " << secondSolutionCentralInstant << "\n\n";
//        std::cout << "thirdSolutionCentralInstant: " << thirdSolutionCentralInstant << "\n\n";


        // Compute partials of central instant w.r.t. link ends position.

        if ( ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) &&
             ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
        {
//            std::cout << "FIRST SOLUTION CENTRAL INSTANT" << "\n\n";
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
//            std::cout << "SECOND SOLUTION CENTRAL INSTANT" << "\n\n";
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
//            std::cout << "THIRD SOLUTION CENTRAL INSTANT" << "\n\n";
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

//        std::cout << "BETA POSITIVE" << "\n\n";

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
        computePartialOfIntermediateVariableTWrtLinkEndState(
                    intermediateVariableQ, intermediateVariableR, intermediateVariableT, partialsIntermediateVariableQwrtFirstTransmitterPosition,
                    partialsIntermediateVariableQwrtSecondTransmitterPosition, partialsIntermediateVariableQwrtReceiverPosition,
                    partialsIntermediateVariableRwrtFirstTransmitterPosition, partialsIntermediateVariableRwrtSecondTransmitterPosition,
                    partialsIntermediateVariableRwrtReceiverPosition, partialsIntermediateVariableTwrtFirstTransmitterPosition,
                    partialsIntermediateVariableTwrtSecondTransmitterPosition, partialsIntermediateVariableTwrtReceiverPosition );

        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtFirstTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtSecondTransmitterPosition;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtReceiverPosition;
        computePartialOfIntermediateVariableSWrtLinkEndState(
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



void MutualApproximationScalingBase::computePartialsOfCentralInstantWrtLinkEndVelocity(
        const std::vector< double > times )
{

    double intermediateVariableQ = ( 3.0 * depressedCubicPolynomialCoefficients_[ 1 ]
            - depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] ) / 9.0;
    double intermediateVariableR = ( 9.0 * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 1 ]
            - 27.0 * depressedCubicPolynomialCoefficients_[ 2 ]
            - 2.0 * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] * depressedCubicPolynomialCoefficients_[ 0 ] ) / 54.0;

    double intermediateVariableBeta = intermediateVariableQ * intermediateVariableQ * intermediateVariableQ
            + intermediateVariableR * intermediateVariableR;

    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity;
    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity;
    Eigen::Matrix< double, 3, 3 > partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity;
    computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianVelocity( partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity,
                                                                              partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity,
                                                                              partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity );

    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtFirstTransmitterVelocity;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtSecondTransmitterVelocity;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableQwrtReceiverVelocity;
    computePartialOfIntermediateVariableQWrtLinkEndState(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity, partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity,
                partialsIntermediateVariableQwrtFirstTransmitterVelocity, partialsIntermediateVariableQwrtSecondTransmitterVelocity,
                partialsIntermediateVariableQwrtReceiverVelocity );

    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtFirstTransmitterVelocity;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtSecondTransmitterVelocity;
    Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableRwrtReceiverVelocity;
    computePartialOfIntermediateVariableRWrtLinkEndState(
                depressedCubicPolynomialCoefficients_, partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity,
                partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity, partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity,
                partialsIntermediateVariableRwrtFirstTransmitterVelocity, partialsIntermediateVariableRwrtSecondTransmitterVelocity,
                partialsIntermediateVariableRwrtReceiverVelocity );

    if ( intermediateVariableBeta < 0 )
    {
//        std::cout << "BETA NEGATITVE" << "\n\n";

        double thetaAngle = computeAngleThetaRealSolutionsCubicEquation( intermediateVariableQ, intermediateVariableR );

        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtFirstTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtSecondTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialOfAngleThetaCubicEquationWrtReceiverVelocity;
        computePartialOfAngleThetaCubicEquationWrtLinkEndState( intermediateVariableQ, intermediateVariableR,
                                                                partialsIntermediateVariableQwrtFirstTransmitterVelocity,
                                                                partialsIntermediateVariableQwrtSecondTransmitterVelocity,
                                                                partialsIntermediateVariableQwrtReceiverVelocity,
                                                                partialsIntermediateVariableRwrtFirstTransmitterVelocity,
                                                                partialsIntermediateVariableRwrtSecondTransmitterVelocity,
                                                                partialsIntermediateVariableRwrtReceiverVelocity,
                                                                partialOfAngleThetaCubicEquationWrtFirstTransmitterVelocity,
                                                                partialOfAngleThetaCubicEquationWrtSecondTransmitterVelocity,
                                                                partialOfAngleThetaCubicEquationWrtReceiverVelocity );

        double firstSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( thetaAngle / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;
        double secondSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;
        double thirdSolutionCentralInstant = 2.0 * std::sqrt( - intermediateVariableQ ) * std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                - depressedCubicPolynomialCoefficients_[ 0 ] / 3.0;

//        std::cout << "estimated central instant: " << times[ 0 ] << "\n\n";
//        std::cout << "firstSolutionCentralInstant: " << firstSolutionCentralInstant << "\n\n";
//        std::cout << "secondSolutionCentralInstant: " << secondSolutionCentralInstant << "\n\n";
//        std::cout << "thirdSolutionCentralInstant: " << thirdSolutionCentralInstant << "\n\n";


        // Compute partials of central instant w.r.t. link ends velocity.

        if ( ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) &&
             ( std::fabs( firstSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
        {
//            std::cout << "FIRST SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterVelocity_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterVelocity_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverVelocity_ =
                    std::cos( ( thetaAngle ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity.block( 0, 0, 1, 3 );
        }
        else if ( ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
                  ( std::fabs( secondSolutionCentralInstant ) <= std::fabs( thirdSolutionCentralInstant ) ) )
        {
//            std::cout << "SECOND SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterVelocity_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterVelocity_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverVelocity_ =
                    std::cos( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 2.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity.block( 0, 0, 1, 3 );
        }
        else if ( ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( firstSolutionCentralInstant ) ) &&
                  ( std::fabs( thirdSolutionCentralInstant ) <= std::fabs( secondSolutionCentralInstant ) ) )
        {
//            std::cout << "THIRD SOLUTION CENTRAL INSTANT" << "\n\n";
            partialsOfCentralInstantWrtFirstTransmitterVelocity_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtFirstTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtFirstTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtSecondTransmitterVelocity_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtSecondTransmitterVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtSecondTransmitterVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 0, 0, 1, 3 );

            partialsOfCentralInstantWrtReceiverVelocity_ =
                    std::cos( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 ) / sqrt( - intermediateVariableQ )
                    * ( - partialsIntermediateVariableQwrtReceiverVelocity )
                    - 2.0 / 3.0 * std::sqrt( - intermediateVariableQ ) * std::sin( ( thetaAngle + 4.0 * mathematical_constants::PI ) / 3.0 )
                    * partialOfAngleThetaCubicEquationWrtReceiverVelocity
                    - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity.block( 0, 0, 1, 3 );
        }

    }

    else
    {

//        std::cout << "BETA POSITIVE" << "\n\n";

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

        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtFirstTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtSecondTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableTwrtReceiverVelocity;
        computePartialOfIntermediateVariableTWrtLinkEndState(
                    intermediateVariableQ, intermediateVariableR, intermediateVariableT, partialsIntermediateVariableQwrtFirstTransmitterVelocity,
                    partialsIntermediateVariableQwrtSecondTransmitterVelocity, partialsIntermediateVariableQwrtReceiverVelocity,
                    partialsIntermediateVariableRwrtFirstTransmitterVelocity, partialsIntermediateVariableRwrtSecondTransmitterVelocity,
                    partialsIntermediateVariableRwrtReceiverVelocity, partialsIntermediateVariableTwrtFirstTransmitterVelocity,
                    partialsIntermediateVariableTwrtSecondTransmitterVelocity, partialsIntermediateVariableTwrtReceiverVelocity );

        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtFirstTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtSecondTransmitterVelocity;
        Eigen::Matrix< double, 1, 3 > partialsIntermediateVariableSwrtReceiverVelocity;
        computePartialOfIntermediateVariableSWrtLinkEndState(
                    intermediateVariableQ, intermediateVariableR, intermediateVariableS, partialsIntermediateVariableQwrtFirstTransmitterVelocity,
                    partialsIntermediateVariableQwrtSecondTransmitterVelocity, partialsIntermediateVariableQwrtReceiverVelocity,
                    partialsIntermediateVariableRwrtFirstTransmitterVelocity, partialsIntermediateVariableRwrtSecondTransmitterVelocity,
                    partialsIntermediateVariableRwrtReceiverVelocity, partialsIntermediateVariableSwrtFirstTransmitterVelocity,
                    partialsIntermediateVariableSwrtSecondTransmitterVelocity, partialsIntermediateVariableSwrtReceiverVelocity );

        // Compute partials of central instant t0 w.r.t. link end velocity.
        partialsOfCentralInstantWrtFirstTransmitterVelocity_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtFirstTransmitterVelocity.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtFirstTransmitterVelocity + partialsIntermediateVariableTwrtFirstTransmitterVelocity;

        partialsOfCentralInstantWrtSecondTransmitterVelocity_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtSecondTransmitterVelocity.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtSecondTransmitterVelocity + partialsIntermediateVariableTwrtSecondTransmitterVelocity;

        partialsOfCentralInstantWrtReceiverVelocity_ =
                - 1.0 / 3.0 * partialOfDepressedPolynomialCoefficientsWrtReceiverVelocity.block( 0, 0, 1, 3 )
                + partialsIntermediateVariableSwrtReceiverVelocity + partialsIntermediateVariableTwrtReceiverVelocity;
    }
}



//! Function to check whether all the required dependent variables are included in the interface.
void MutualApproximationScalingBase::checkRequiredDependentVariablesFromInterface( const observation_models::LinkEnds linkEnds )
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
void MutualApproximationScalingBase::retrieveRelativeAccelerationsAndAssociatedPartialsFromDependentVariables(
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
        partialAccelerationFirstTransmitterWrtReceiverVelocity_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtReceiverPosition.segment( i * 6 + 3, 3 ).transpose( );
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
        partialAccelerationFirstTransmitterWrtTransmitterVelocity_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtTransmitterPosition.segment( i * 6 + 3, 3 ).transpose( );
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
        partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity_.block( i, 0, 1, 3 ) = partialAccelerationFirstTransmitterWrtOtherTransmitterPosition.segment( i * 6 + 3, 3 ).transpose( );
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
        partialAccelerationSecondTransmitterWrtReceiverVelocity_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtReceiverPosition.segment( i * 6 + 3, 3 ).transpose( );
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
        partialAccelerationSecondTransmitterWrtTransmitterVelocity_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtTransmitterPosition.segment( i * 6 + 3, 3 ).transpose( );
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
        partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity_.block( i, 0, 1, 3 ) = partialAccelerationSecondTransmitterWrtOtherTransmitterPosition.segment( i * 6 + 3, 3 ).transpose( );
    }
}

//! Update the scaling object to the current times and states
void MutualApproximationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                         const std::vector< double >& times,
                                         const observation_models::LinkEndType fixedLinkEnd,
                                         const observation_models::LinkEnds linkEnds,
                                         const Eigen::VectorXd currentObservation )
{
//    std::cout << "BEGINNING UPDATE FUNCTION IN MUTUAL APPROXIMATION SCALING" << "\n\n";

//    std::cout << "first transmitter time: " << times[ 0 ] << "\n\n";
//    std::cout << "second transmitter time: " << times[ 1 ] << "\n\n";
//    std::cout << "receiver time: " << times[ 2 ] << "\n\n";

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

    instrumentalFrameRelativeAcceleration_ = computeRelativeAccelerationInInstrumentalFrame( );
//    std::cout << "instrumental frame relative acceleration: " << instrumentalFrameRelativeAcceleration_.transpose( ) << "\n\n";



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian position.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPositionV2( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPositionV2( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian positions.
//    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
//                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
//                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
//                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_ );
//    std::cout << "partials of relative acceleration w.r.t. first transmitter position (v1): " << partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partials of relative acceleration w.r.t. second transmitter position (v1): " << partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partials of relative acceleration w.r.t. receiver position (v1): " << partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ << "\n\n";

    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPositionV2(
                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_ );
//    std::cout << "partials of relative acceleration w.r.t. first transmitter position (v2): " << partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partials of relative acceleration w.r.t. second transmitter position (v2): " << partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partials of relative acceleration w.r.t. receiver position (v2): " << partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ << "\n\n";


    // Compute coefficients of the cubic polynomial for the central instant t0.
    cubicPolynomialCoefficients_ = computeCubicPolynomialCoefficients( );

    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
    depressedCubicPolynomialCoefficients_ = computeDepressedCubicPolynomialCoefficients( );

    // Compute partials of central instant w.r.t. link ends positions.
    computePartialsOfCentralInstantWrtLinkEndPosition( times );



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndVelocity( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndVelocity( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndVelocity(
                partialAccelerationFirstTransmitterWrtReceiverVelocity_,
                partialAccelerationFirstTransmitterWrtTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtReceiverVelocity_,
                partialAccelerationSecondTransmitterWrtTransmitterVelocity_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity_ );

//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition: " << partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition: " << partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition: " << partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_ << "\n\n";

//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition: " << partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition: " << partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition: " << partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_ << "\n\n";

//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition: " << partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ << "\n\n";


//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterVelocity: " << partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterVelocity: " << partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativePositionWrtReceiverVelocity: " << partialsOfInstrumentalFrameRelativePositionWrtReceiverVelocity_ << "\n\n";

//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity: " << partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity: " << partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity: " << partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity_ << "\n\n";

//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterVelocity: " << partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterVelocity: " << partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterVelocity_ << "\n\n";
//    std::cout << "partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverVelocity: " << partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverVelocity_ << "\n\n";


    // Compute partials of central instant w.r.t. link ends velocities.
    computePartialsOfCentralInstantWrtLinkEndVelocity( times );


    currentLinkEndType_ = fixedLinkEnd;

}



void MutualApproximationWithImpactParameterScaling::computePartialsOfApparentDistanceWrtLinkEndPosition( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;
        }


        // Compute partials of apparent distance w.r.t. position of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  ( 1.0 / apparentDistance ) *
                ( instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition.block( 1, 0, 1, 3 ) );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfApparentDistanceWrtFirstTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfApparentDistanceWrtSecondTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfApparentDistanceWrtReceiverPosition_ = partials;
        }

    }
}


void MutualApproximationWithImpactParameterScaling::computePartialsOfImpactParameterWrtLinkEndPosition( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    double timeDerivativeOfApparentDistance = ( 1.0 / apparentDistance )
            * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfApparentDistanceWrtLinkEndPosition;
        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtFirstTransmitterPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtSecondTransmitterPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtReceiverPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtReceiverPosition_;
        }


        // Compute partials of impact parameter w.r.t. position of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  partialsOfApparentDistanceWrtLinkEndPosition
                + timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndPosition;


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfImpactParameterWrtFirstTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfImpactParameterWrtSecondTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfImpactParameterWrtReceiverPosition_ = partials;
        }

    }
}



void MutualApproximationWithImpactParameterScaling::computePartialsOfImpactParameterWrtLinkEndVelocity( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    double timeDerivativeOfApparentDistance = ( 1.0 / apparentDistance )
            * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtLinkEndVelocity;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtFirstTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtSecondTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtReceiverVelocity_;
        }


        // Compute partials of impact parameter w.r.t. velocity of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndVelocity;


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfImpactParameterWrtFirstTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfImpactParameterWrtSecondTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver velocity
        {
            partialsOfImpactParameterWrtReceiverVelocity_ = partials;
        }

    }
}



//! Update the scaling object to the current times and states
void MutualApproximationWithImpactParameterScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                            const std::vector< double >& times,
                                                            const observation_models::LinkEndType fixedLinkEnd,
                                                            const observation_models::LinkEnds linkEnds,
                                                            const Eigen::VectorXd currentObservation )
{
//    std::cout << "BEGINNING UPDATE FUNCTION IN MUTUAL APPROXIMATION WITH IMPACT PARAMETER SCALING" << "\n\n";

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

    instrumentalFrameRelativeAcceleration_ = computeRelativeAccelerationInInstrumentalFrame( );
//    std::cout << "instrumental frame relative acceleration: " << instrumentalFrameRelativeAcceleration_.transpose( ) << "\n\n";



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian position.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPositionV2( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPositionV2( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian positions.
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPositionV2(
                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_ );

    // Compute coefficients of the cubic polynomial for the central instant t0.
    cubicPolynomialCoefficients_ = computeCubicPolynomialCoefficients( );

    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
    depressedCubicPolynomialCoefficients_ = computeDepressedCubicPolynomialCoefficients( );

    // Compute partials of central instant w.r.t. link ends positions.
    computePartialsOfCentralInstantWrtLinkEndPosition( times );


    // Compute partials of apparent distance w.r.t. link ends positions.
    computePartialsOfApparentDistanceWrtLinkEndPosition( );

    // Compute partials of impact parameter w.r.t. link ends position.
    computePartialsOfImpactParameterWrtLinkEndPosition( );


    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndVelocity( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndVelocity( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndVelocity(
                partialAccelerationFirstTransmitterWrtReceiverVelocity_,
                partialAccelerationFirstTransmitterWrtTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtReceiverVelocity_,
                partialAccelerationSecondTransmitterWrtTransmitterVelocity_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity_ );

    // Compute partials of central instant w.r.t. link ends velocities.
    computePartialsOfCentralInstantWrtLinkEndVelocity( times );

    // Compute partials of impact parameter w.r.t. link ends velocity.
    computePartialsOfImpactParameterWrtLinkEndVelocity( );



    currentLinkEndType_ = fixedLinkEnd;

}



void ModifiedMutualApproximationScaling::computePartialsOfModifiedObservableWrtLinkEndPosition( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 2, 3 > partialsOfRelativePositionWrtLinkEndPosition;
        Eigen::Matrix< double, 2, 3 > partialsOfRelativeVelocityWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;
            partialsOfRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;
            partialsOfRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;
            partialsOfRelativeVelocityWrtLinkEndPosition = partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_;
        }


        // Compute partials of impact parameter w.r.t. position of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  ( 1.0 / apparentDistance )
                * ( instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativeVelocityWrtLinkEndPosition.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 1, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativeVelocityWrtLinkEndPosition.block( 1, 0, 1, 3 ) )
                - ( 1.0 / ( apparentDistance * apparentDistance * apparentDistance ) )
                * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
                + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] )
                * ( instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 1, 0, 1, 3 ) );
//                ( instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 0, 0, 1, 3 )
//                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativeVelocityWrtLinkEndPosition.block( 0, 0, 1, 3 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfRelativePositionWrtLinkEndPosition.block( 1, 0, 1, 3 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativeVelocityWrtLinkEndPosition.block( 1, 0, 1, 3 ) );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfModifiedObservableWrtFirstTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfModifiedObservableWrtSecondTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfModifiedObservableWrtReceiverPosition_ = partials;
        }

    }
}



void ModifiedMutualApproximationScaling::computePartialsOfModifiedObservableWrtLinkEndVelocity( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 2, 3 > partialsOfRelativePositionWrtLinkEndVelocity;
        Eigen::Matrix< double, 2, 3 > partialsOfRelativeVelocityWrtLinkEndVelocity;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterVelocity_;
            partialsOfRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterVelocity_;
            partialsOfRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver velocity
        {
            partialsOfRelativePositionWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativePositionWrtReceiverVelocity_;
            partialsOfRelativeVelocityWrtLinkEndVelocity = partialsOfInstrumentalFrameRelativeVelocityWrtReceiverVelocity_;
        }

        // Compute partials of impact parameter w.r.t. velocity of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  ( 1.0 / apparentDistance )
                * ( instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativeVelocityWrtLinkEndVelocity.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 1, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativeVelocityWrtLinkEndVelocity.block( 1, 0, 1, 3 ) )
                - ( 1.0 / ( apparentDistance * apparentDistance * apparentDistance ) )
                * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
                + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] )
                * ( instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 1, 0, 1, 3 ) );
//                ( instrumentalFrameRelativeVelocity_[ 0 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 0, 0, 1, 3 )
//                + instrumentalFrameRelativePosition_[ 0 ] * partialsOfRelativeVelocityWrtLinkEndVelocity.block( 0, 0, 1, 3 )
//                + instrumentalFrameRelativeVelocity_[ 1 ] * partialsOfRelativePositionWrtLinkEndVelocity.block( 1, 0, 1, 3 )
//                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfRelativeVelocityWrtLinkEndVelocity.block( 1, 0, 1, 3 ) );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfModifiedObservableWrtFirstTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfModifiedObservableWrtSecondTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver velocity
        {
            partialsOfModifiedObservableWrtReceiverVelocity_ = partials;
        }

    }
}


//! Update the scaling object to the current times and states
void ModifiedMutualApproximationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                 const std::vector< double >& times,
                                                 const observation_models::LinkEndType fixedLinkEnd,
                                                 const observation_models::LinkEnds linkEnds,
                                                 const Eigen::VectorXd currentObservation )
{
//    std::cout << "BEGINNING UPDATE FUNCTION IN MODIFIED MUTUAL APPROXIMATION SCALING" << "\n\n";

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


//    // Compute second partials of right ascension and declination of the two transmitters wrt time.
//    secondPartialOfRightAscensionFirstTransmitterWrtTime_ = computeSecondPartialRightAscensionWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationFirstTransmitterWrtReceiver_ );

//    secondPartialOfRightAscensionSecondTransmitterWrtTime_ = computeSecondPartialRightAscensionWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationSecondTransmitterWrtReceiver_ );

//    secondPartialOfDeclinationFirstTransmitterWrtTime_ = computeSecondPartialDeclinationWrtTime(
//                ( firstTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( firstTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationFirstTransmitterWrtReceiver_ );

//    secondPartialOfDeclinationSecondTransmitterWrtTime_ = computeSecondPartialDeclinationWrtTime(
//                ( secondTransmitterState_ - receiverState_ ).segment( 0, 3 ),
//                ( secondTransmitterState_ - receiverState_ ).segment( 3, 3 ),
//                cartesianAccelerationSecondTransmitterWrtReceiver_ );


    // Compute relative position, velocity and acceleration between two transmitters in the instrumental frame of
    // the receiver.
    instrumentalFrameRelativePosition_ = computeRelativePositionInInstrumentalFrame( );
//    std::cout << "instrumental frame relative position: " << instrumentalFrameRelativePosition_.transpose( ) << "\n\n";

    instrumentalFrameRelativeVelocity_ = computeRelativeVelocityInInstrumentalFrame( );
//    std::cout << "instrumental frame relative velocity: " << instrumentalFrameRelativeVelocity_.transpose( ) << "\n\n";

//    instrumentalFrameRelativeAcceleration_ = computeRelativeAccelerationInInstrumentalFrame( );
////    std::cout << "instrumental frame relative acceleration: " << instrumentalFrameRelativeAcceleration_.transpose( ) << "\n\n";


//    // Compute apparent distance.
//    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
//            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian position.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPositionV2( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPositionV2( );

    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndVelocity( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndVelocity( );

//    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian positions.
//    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
//                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
//                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
//                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_,
//                partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_, partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_,
//                partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_ );

//    // Compute coefficients of the cubic polynomial for the central instant t0.
//    cubicPolynomialCoefficients_ = computeCubicPolynomialCoefficients( );

//    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
//    depressedCubicPolynomialCoefficients_ = computeDepressedCubicPolynomialCoefficients( );

//    // Compute partials of central instant w.r.t. link ends positions.
//    computePartialsOfCentralInstantWrtLinkEndPosition( times );

    // Compute partial of the modified mutual approximation observable w.r.t. link end cartesian positions.
    computePartialsOfModifiedObservableWrtLinkEndPosition( );

    // Compute partial of the modified mutual approximation observable w.r.t. link end cartesian velocities.
    computePartialsOfModifiedObservableWrtLinkEndVelocity( );


    currentLinkEndType_ = fixedLinkEnd;

}




void ImpactParameterMutualApproxScaling::computePartialsOfApparentDistanceWrtLinkEndPosition( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition = partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;
        }


        // Compute partials of apparent distance w.r.t. position of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  ( 1.0 / apparentDistance ) *
                ( instrumentalFrameRelativePosition_[ 0 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition.block( 0, 0, 1, 3 )
                + instrumentalFrameRelativePosition_[ 1 ] * partialsOfInstrumentalFrameRelativePositionWrtLinkEndPosition.block( 1, 0, 1, 3 ) );


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfApparentDistanceWrtFirstTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfApparentDistanceWrtSecondTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfApparentDistanceWrtReceiverPosition_ = partials;
        }

    }
}


void ImpactParameterMutualApproxScaling::computePartialsOfImpactParameterWrtLinkEndPosition( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    double timeDerivativeOfApparentDistance = ( 1.0 / apparentDistance )
            * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfApparentDistanceWrtLinkEndPosition;
        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtLinkEndPosition;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtFirstTransmitterPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtFirstTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtSecondTransmitterPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtSecondTransmitterPosition_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver position
        {
            partialsOfApparentDistanceWrtLinkEndPosition = partialsOfApparentDistanceWrtReceiverPosition_;
            partialsOfCentralInstantWrtLinkEndPosition = partialsOfCentralInstantWrtReceiverPosition_;
        }


        // Compute partials of impact parameter w.r.t. position of current link end .
        Eigen::Matrix< double, 1, 3 > partials =  partialsOfApparentDistanceWrtLinkEndPosition
                + timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndPosition;

        std::cout << "partials apparent distance w.r.t. position: " << partialsOfApparentDistanceWrtLinkEndPosition << "\n\n";
        std::cout << "partials contribution central instant: " << timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndPosition << "\n\n";


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter position
        {
            partialsOfImpactParameterWrtFirstTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter position
        {
            partialsOfImpactParameterWrtSecondTransmitterPosition_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver position
        {
            partialsOfImpactParameterWrtReceiverPosition_ = partials;
        }

    }
}



void ImpactParameterMutualApproxScaling::computePartialsOfImpactParameterWrtLinkEndVelocity( )
{
    double apparentDistance = std::sqrt( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativePosition_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativePosition_[ 1 ] );

    double timeDerivativeOfApparentDistance = ( 1.0 / apparentDistance )
            * ( instrumentalFrameRelativePosition_[ 0 ] * instrumentalFrameRelativeVelocity_[ 0 ]
            + instrumentalFrameRelativePosition_[ 1 ] * instrumentalFrameRelativeVelocity_[ 1 ] );


    // Successively compute partials w.r.t. state of first transmitter, second transmitter and receiver.
    for ( unsigned int currentLinkEndIndex = 0 ; currentLinkEndIndex < 3 ; currentLinkEndIndex++ )
    {

        Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtLinkEndVelocity;
        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtFirstTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtSecondTransmitterVelocity_;
        }
        else if ( currentLinkEndIndex == 2 ) // If partials w.r.t. receiver velocity
        {
            partialsOfCentralInstantWrtLinkEndVelocity = partialsOfCentralInstantWrtReceiverVelocity_;
        }


        // Compute partials of impact parameter w.r.t. velocity of current link end .
        Eigen::Matrix< double, 1, 3 > partials = timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndVelocity;

        std::cout << "partials apparent distance w.r.t. velocity: " << timeDerivativeOfApparentDistance * partialsOfCentralInstantWrtLinkEndVelocity << "\n\n";


        // Assign partials to the right link leg.

        if ( currentLinkEndIndex == 0 ) // If partials w.r.t. first transmitter velocity
        {
            partialsOfImpactParameterWrtFirstTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 1 ) // If partials w.r.t. second transmitter velocity
        {
            partialsOfImpactParameterWrtSecondTransmitterVelocity_ = partials;
        }
        else if ( currentLinkEndIndex == 2 )  // If partials w.r.t. receiver velocity
        {
            partialsOfImpactParameterWrtReceiverVelocity_ = partials;
        }

    }
}



//! Update the scaling object to the current times and states
void ImpactParameterMutualApproxScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                 const std::vector< double >& times,
                                                 const observation_models::LinkEndType fixedLinkEnd,
                                                 const observation_models::LinkEnds linkEnds,
                                                 const Eigen::VectorXd currentObservation )
{
//    std::cout << "BEGINNING UPDATE FUNCTION IN MUTUAL APPROXIMATION WITH IMPACT PARAMETER SCALING" << "\n\n";

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

    instrumentalFrameRelativeAcceleration_ = computeRelativeAccelerationInInstrumentalFrame( );
//    std::cout << "instrumental frame relative acceleration: " << instrumentalFrameRelativeAcceleration_.transpose( ) << "\n\n";



    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian position.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPositionV2( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPositionV2( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian positions.
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPositionV2(
                cartesianAccelerationFirstTransmitterWrtReceiver_, cartesianAccelerationSecondTransmitterWrtReceiver_, partialAccelerationFirstTransmitterWrtReceiverPosition_,
                partialAccelerationFirstTransmitterWrtTransmitterPosition_, partialAccelerationSecondTransmitterWrtReceiverPosition_, partialAccelerationSecondTransmitterWrtTransmitterPosition_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_, partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_ );

    // Compute coefficients of the cubic polynomial for the central instant t0.
    cubicPolynomialCoefficients_ = computeCubicPolynomialCoefficients( );

    // Compute coefficients of the depressed cubic polynomial for the central instant t0.
    depressedCubicPolynomialCoefficients_ = computeDepressedCubicPolynomialCoefficients( );

    // Compute partials of central instant w.r.t. link ends positions.
    computePartialsOfCentralInstantWrtLinkEndPosition( times );


    // Compute partials of apparent distance w.r.t. link ends positions.
    computePartialsOfApparentDistanceWrtLinkEndPosition( );

    // Compute partials of impact parameter w.r.t. link ends position.
    computePartialsOfImpactParameterWrtLinkEndPosition( );


    // Compute partials of relative position and velocity between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndVelocity( );
    computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndVelocity( );

    // Compute partials of relative acceleration between two transmitters in the instrumental frame wrt link end cartesian velocity.
    computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndVelocity(
                partialAccelerationFirstTransmitterWrtReceiverVelocity_,
                partialAccelerationFirstTransmitterWrtTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtReceiverVelocity_,
                partialAccelerationSecondTransmitterWrtTransmitterVelocity_,
                partialAccelerationFirstTransmitterWrtOtherTransmitterVelocity_,
                partialAccelerationSecondTransmitterWrtOtherTransmitterVelocity_ );

    // Compute partials of central instant w.r.t. link ends velocities.
    computePartialsOfCentralInstantWrtLinkEndVelocity( times );

    // Compute partials of impact parameter w.r.t. link ends velocity.
    computePartialsOfImpactParameterWrtLinkEndVelocity( );


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
        throw std::runtime_error( "Error mutual approximation partial and scaling are inconsistent." );
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


        // Retrieve scaling factor.
        Eigen::Matrix< double, 1, 3 > scalingFactor = Eigen::Matrix< double, 1, 3 >::Zero( );
        Eigen::Matrix< double, 1, 3 > velocityScalingFactor = Eigen::Matrix< double, 1, 3 >::Zero( );
        if ( std::dynamic_pointer_cast< MutualApproximationScaling >( mutualApproximationScaler_ ) != nullptr )
        {
            scalingFactor = std::dynamic_pointer_cast< MutualApproximationScaling >( mutualApproximationScaler_ )
                    ->getScalingFactor( positionPartialIterator_->first );
            velocityScalingFactor = std::dynamic_pointer_cast< MutualApproximationScaling >( mutualApproximationScaler_ )
                    ->getVelocityScalingFactor( positionPartialIterator_->first );
        }
        else
        {
            scalingFactor = std::dynamic_pointer_cast< ModifiedMutualApproximationScaling >( mutualApproximationScaler_ )
                    ->getScalingFactor( positionPartialIterator_->first );
            velocityScalingFactor = std::dynamic_pointer_cast< ModifiedMutualApproximationScaling >( mutualApproximationScaler_ )
                    ->getVelocityScalingFactor( positionPartialIterator_->first );
        }

        // Scale position partials
        returnPartial.push_back( std::make_pair(
                                     scalingFactor * ( positionPartialIterator_->second->calculatePartialOfPosition( currentState_ , currentTime_ ) )
                                     + velocityScalingFactor * ( positionPartialIterator_->second->calculatePartialOfVelocity( currentState_ , currentTime_ ) ),
                                     currentTime_ ) );

//        std::cout << "TEST returnPartial state partial: " <<
//                     scalingFactor * ( positionPartialIterator_->second->calculatePartialOfPosition( currentState_ , currentTime_ ) )
//                     + velocityScalingFactor * ( positionPartialIterator_->second->calculatePartialOfVelocity( currentState_ , currentTime_ ) ) << "\n\n";


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



//! Function to calculate the observation partial(s) at required time and state
MutualApproximationWithImpactParameterPartial::MutualApproximationWithImpactParameterPartialReturnType MutualApproximationWithImpactParameterPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector2d& currentObservation )
{
    if( linkEndOfFixedTime != mutualApproximationScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error mutual approximation with impact parameter partial and scaling are inconsistent." );
    }

    MutualApproximationWithImpactParameterPartialReturnType returnPartial;

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
                        mutualApproximationScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) )
                        + mutualApproximationScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfVelocity(
                              currentState_ , currentTime_ ) ), currentTime_ ) );

//        std::cout << "TEST returnPartial state partial: " <<
//                     mutualApproximationScaler_->getScalingFactor( positionPartialIterator_->first ) *
//                     ( positionPartialIterator_->second->calculatePartialOfPosition(
//                           currentState_ , currentTime_ ) )
//                     + mutualApproximationScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
//                     ( positionPartialIterator_->second->calculatePartialOfVelocity(
//                           currentState_ , currentTime_ ) ) << "\n\n";

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


//! Function to calculate the observation partial(s) at required time and state
ImpactParameterMutualApproxPartial::ImpactParameterMutualApproxPartialReturnType ImpactParameterMutualApproxPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    if( linkEndOfFixedTime != impactParameterScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error impact parameter (mutual approximation) partial and scaling are inconsistent." );
    }

    ImpactParameterMutualApproxPartialReturnType returnPartial;

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
                        impactParameterScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) )
                        + impactParameterScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfVelocity(
                              currentState_ , currentTime_ ) ), currentTime_ ) );

//        std::cout << "TEST returnPartial state partial: " <<
//                     impactParameterScaler_->getScalingFactor( positionPartialIterator_->first ) *
//                     ( positionPartialIterator_->second->calculatePartialOfPosition(
//                           currentState_ , currentTime_ ) )
//                     + impactParameterScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
//                     ( positionPartialIterator_->second->calculatePartialOfVelocity(
//                           currentState_ , currentTime_ ) ) << "\n\n";

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
