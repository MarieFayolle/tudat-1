/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/inertialTorquePartial.h"

namespace tudat
{

namespace acceleration_partials
{

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
InertialTorquePartial::getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >  partialFunction = std::make_pair(
                boost::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case mean_moment_of_inertia:
        {
            partialFunction = std::make_pair(
                        boost::bind( &InertialTorquePartial::wrtMeanMomentOfInertia, this, _1 ), 1 );
            break;
        }
        default:
            break;
        }
    }
    return partialFunction;

}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
/*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >  partialFunction = std::make_pair(
                boost::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case spherical_harmonics_cosine_coefficient_block:
        {
            // Cast parameter object to required type.
            boost::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                    boost::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

            int c20Index, c21Index, c22Index;
            coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );

            if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
            {
                if( getInertiaTensorNormalizationFactor_.empty( ) )
                {
                    throw std::runtime_error( "Error when getting partial of 2nd degree grac torque w.r.t. cosine sh parameters, inertia tensor normalization function not found." );
                }
                partialFunction = std::make_pair(
                            boost::bind( &InertialTorquePartial::
                                         wrtCosineSphericalHarmonicCoefficientsOfCentralBody, this,
                                         _1, c20Index, c21Index, c22Index ), coefficientsParameter->getParameterSize( ) );
            }

            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            // Cast parameter object to required type.

            boost::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                    boost::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

            int s21Index, s22Index;
            coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );

            if( s21Index >= 0 || s22Index >= 0 )
            {
                if( getInertiaTensorNormalizationFactor_.empty( ) )
                {
                    throw std::runtime_error( "Error when getting partial of 2nd degree grac torque w.r.t. sine sh parameters, inertia tensor normalization function not found." );
                }
                partialFunction = std::make_pair(
                            boost::bind( &InertialTorquePartial::
                                         wrtSineSphericalHarmonicCoefficientsOfCentralBody, this,
                                         _1, s21Index, s22Index ), coefficientsParameter->getParameterSize( ) );
            }


            break;
        }
        default:
            break;
        }
    }
    return partialFunction;
}

void InertialTorquePartial::wrtMeanMomentOfInertia(
        Eigen::MatrixXd& momentOfInertiaPartial )
{
    momentOfInertiaPartial .block( 0, 1, 3, 1 ) =
            UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_MEAN_MOMENT *
            currentInverseInertiaTensor_* currentTotalTorque_;
}

void InertialTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody(
        Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
        const int c20Index, const int c21Index, const int c22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( c20Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c20Index, 3, 1 ) =
        UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 *
                currentInverseInertiaTensor_* currentTotalTorque_;
    }

    if( c21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c21Index, 3, 1 ) =
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 *
                currentInverseInertiaTensor_* currentTotalTorque_;
    }

    if( c22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c22Index, 3, 1 ) =
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 *
                currentInverseInertiaTensor_* currentTotalTorque_;
    }
}

void InertialTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody(
        Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
        const int s21Index, const int s22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( s21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, s21Index, 3, 1 ) =
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 *
                currentInverseInertiaTensor_* currentTotalTorque_;
    }

    if( s22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, s22Index, 3, 1 ) =
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 *
                currentInverseInertiaTensor_* currentTotalTorque_;
    }
}


} // namespace acceleration_partials

} // namespace tudat
