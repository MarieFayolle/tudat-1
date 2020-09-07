/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALAPPROXIMATIONPARTIAL_H
#define TUDAT_MUTUALAPPROXIMATIONPARTIAL_H

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/angularPositionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/apparentDistancePartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"

namespace tudat
{

namespace observation_partials
{

double computePartialOfRightAscensionWrtTime( Eigen::Vector3d cartesianPositionVector,
                                              Eigen::Vector3d cartesianVelocityVector );

double computePartialOfDeclinationWrtTime( Eigen::Vector3d cartesianPositionVector,
                                           Eigen::Vector3d cartesianVelocityVector );

double computeSecondPartialRightAscensionWrtTime( Eigen::Vector3d cartesianPosition,
                                                  Eigen::Vector3d cartesianVelocity,
                                                  Eigen::Vector3d cartesianAcceleration );

double computeSecondPartialDeclinationWrtTime( Eigen::Vector3d cartesianPosition,
                                               Eigen::Vector3d cartesianVelocity,
                                               Eigen::Vector3d cartesianAcceleration );

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector );

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > computePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d cartesianPositionVector );

//! Function to compute the partial w.r.t. position of observer or observed object of the first time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector );


//! Function to compute the partial w.r.t. position of observer or observed object of the first time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfFirstTimeDerivativeDeclinationWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector );

//! Function to compute the partial w.r.t. position of observer or observed object of the second time derivative of the right ascension.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
        const Eigen::Vector3d& cartesianAccelerationVector,
        const Eigen::Matrix3d& partialCartesianAccelerationWrtPosition,
        const bool wrtTransmitterPosition,
        const bool wrtOtherLinkEnd = false );

//! Function to compute the partial w.r.t. position of observer or observed object of the second time derivative of the declination.
Eigen::Matrix< double, 1, 3 > computePartialOfSecondTimeDerivativeDeclinationWrtLinkEndPosition(
        const Eigen::Vector3d& cartesianPositionVector,
        const Eigen::Vector3d& cartesianVelocityVector,
        const Eigen::Vector3d& cartesianAccelerationVector,
        const Eigen::Matrix3d& partialCartesianAccelerationWrtPosition,
        const bool wrtTransmitterPosition,
        const bool wrtOtherLinkEnd = false );

/// TO BE MOVED TO BASIC_MATHEMATICS FRAMEWORK
//! Compute the angle theta used to derive the real solutions of the depressed cubic equation.
double computeAngleThetaRealSolutionsCubicEquation( double intermediateQ,
                                                    double intermediateR );


//! Derived class for scaling three-dimensional position partial to mutual approximation observable partial
class MutualApproximationScalingBase: public PositionPartialScaling
{
public:

    //! Constructor
    MutualApproximationScalingBase(
            const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface ) :
    dependentVariablesInterface_( dependentVariablesInterface )
    {
        dependentVariablesInterfaceIdsAndIndices_ = dependentVariablesInterface_->getDependentVariablesIdsAndIndices( );

        for ( std::map< std::string, int >::iterator itr = dependentVariablesInterfaceIdsAndIndices_.begin( ) ;
              itr != dependentVariablesInterfaceIdsAndIndices_.end( ) ; itr++ )
        {
            std::cout << "dependent variables from interface ID: " << itr->first << "\n\n";
        }
//        referenceAngularPositionScaler_ = std::make_shared< AngularPositionScaling >( );
    }

    //! Destructor
    ~MutualApproximationScalingBase( ){ }

    //! Update the scaling object to the current times and states (pure virtual).
        /*!
         *  Update the scaling object to the current times and states (pure virtual).
         *  \param linkEndStates List of states at each link end during observation.
         *  \param times List of times at each link end during observation.
         *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
         *  is kept constant when computing observable.
         *  \param currentObservation Value of observation for which partial scaling is to be computed
         */
    virtual void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                             const std::vector< double >& times,
                             const observation_models::LinkEndType fixedLinkEnd,
                             const observation_models::LinkEnds linkEnds,
                             const Eigen::VectorXd currentObservation ) = 0;

//    //! Function to retrieve the scaling factor for specific link end
//    /*!
//     * Function to retrieve the scaling factor for specific link end
//     * \param linkEndType Link end for which scaling factor is to be returned
//     * \return Position partial scaling factor at current link end
//     */
//    Eigen::Matrix< double, 1, 3 > getScalingFactor(
//            const observation_models::LinkEndType linkEndType )
//    {
//        if ( linkEndType == observation_models::transmitter )
//        {
//            return partialsOfCentralInstantWrtFirstTransmitterPosition_;
//        }
//        else if ( linkEndType == observation_models::transmitter2 )
//        {
//            return partialsOfCentralInstantWrtSecondTransmitterPosition_;
//        }
//        else if ( linkEndType == observation_models::receiver )
//        {
//            return partialsOfCentralInstantWrtReceiverPosition_;
//        }
//        else
//        {
//            throw std::runtime_error( "Error when retrieving the apparent distance scaling factor, the link end type"
//                                      "is different from transmitter, transmitter2 or receiver" );
//        }
//    }

//    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
//    /*!
//     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
//     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
//     */
//    double getLightTimePartialScalingFactor( )
//    {
////        XlightTimeCorrectionScalingFactor_ = XlightTimeCorrectionScalingFactors_.first
////                + XlightTimeCorrectionScalingFactors_.second;
////        YlightTimeCorrectionScalingFactor_ = YlightTimeCorrectionScalingFactors_.first
////                + YlightTimeCorrectionScalingFactors_.second;

//        return 0.0; //1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
////                * ( XcoordinateReceiverFrame_ * XlightTimeCorrectionScalingFactor_
////                    + YcoordinateReceiverFrame_ * YlightTimeCorrectionScalingFactor_ );
//    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

    //! Function to retrieve the relative position between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d getInstrumentalFrameRelativePosition( )
    {
        return instrumentalFrameRelativePosition_;
    }

    //! Function to retrieve the relative velocity between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d getInstrumentalFrameRelativeVelocity( )
    {
        return instrumentalFrameRelativeVelocity_;
    }

    //! Function to retrieve the relative acceleration between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d getInstrumentalFrameRelativeAcceleration( )
    {
        return instrumentalFrameRelativeAcceleration_;
    }

    //! Function to retrieve the partial of the relative position of the two transmitters in the instrumental frame w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativePositionInInstrumentalFrameWrtFirstTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative position of the two transmitters in the instrumental frame w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativePositionInInstrumentalFrameWrtSecondTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative position of the two transmitters in the instrumental frame w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativePositionInInstrumentalFrameWrtReceiverPosition( )
    {
        return partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;
    }

    //! Function to retrieve the partial of the relative velocity of the two transmitters in the instrumental frame w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativeVelocityInInstrumentalFrameWrtFirstTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative velocity of the two transmitters in the instrumental frame w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativeVelocityInInstrumentalFrameWrtSecondTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative velocity of the two transmitters in the instrumental frame w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfRelativeVelocityInInstrumentalFrameWrtReceiverPosition( )
    {
        return partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_;
    }

    //! Function to retrieve the partial of the relative acceleration of the two transmitters in the instrumental w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative acceleration of the two transmitters in the instrumental w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition( )
    {
        return partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_;
    }

    //! Function to retrieve the partial of the relative acceleration of the two transmitters in the instrumental w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > getPartialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition( )
    {
        return partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_;
    }

protected:

    //! Function to check whether all the required dependent variables are included in the interface.
    void checkRequiredDependentVariablesFromInterface( const observation_models::LinkEnds linkEnds );

    //! Function to retrieve the relative accelerations of the two transmitters, and the associated acceleration partials
    //!  from the dependent variables interface.
    void retrieveRelativeAccelerationsAndAssociatedPartialsFromDependentVariables(
            const std::vector< double >& times,
            const observation_models::LinkEnds linkEnds );

    Eigen::Vector2d computeRelativePositionInInstrumentalFrame( );

    Eigen::Vector2d computeRelativeVelocityInInstrumentalFrame( );

    Eigen::Vector2d computeRelativeAccelerationInInstrumentalFrame( );

    void computePartialOfRelativePositionInInstrumentalFrameWrtLinkEndPosition( );

    void computePartialOfRelativeVelocityInInstrumentalFrameWrtLinkEndPosition( );

    void computePartialOfRelativeAccelerationInInstrumentalFrameWrtLinkEndPosition(
            Eigen::Vector3d relativeAccelerationFirstTransmitterWrtReceiver,
            Eigen::Vector3d relativeAccelerationSecondTransmitterWrtReceiver,
            Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverPosition,
            Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterPosition,
            Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverPosition,
            Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterPosition,
            Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterPosition,
            Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterPosition,
            Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition,
            Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition,
            Eigen::Matrix< double, 2, 3 >& partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition );


    Eigen::Vector4d computeCubicPolynomialCoefficients( );

    void computePartialOfCubicPolynomialCoefficientsWrtCartesianPosition(
            Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
            Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
            Eigen::Matrix< double, 4, 3 >& partialOfCubicPolynomialCoefficientsWrtReceiverPosition );

    Eigen::Vector3d computeDepressedCubicPolynomialCoefficients( );

    void computePartialOfDepressedCubicPolynomialCoefficientsWrtCartesianPosition(
            Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtFirstTransmitterPosition,
            Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtSecondTransmitterPosition,
            Eigen::Matrix< double, 3, 3 >& partialOfDepressedCubicPolynomialCoefficientsWrtReceiverPosition );


    void computePartialsOfCentralInstantWrtLinkEndPosition( const std::vector< double > times );


//private:

    //! Dependent variable interface object, to be used to retrieve total acceleration and total acceleration
    //! partial w.r.t. link end position
    std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface_;

    //! Dependent variables IDs and indices for the dependent variables interface.
    std::map< std::string, int > dependentVariablesInterfaceIdsAndIndices_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;


    //! Computed scaling factor for first transmitter (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorFirstTransmitter_;

    //! Computed scaling factor for second transmitter (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorSecondTransmitter_;


    //!
    std::shared_ptr< AngularPositionScaling > referenceAngularPositionScaler_;

    //! Computed angular position scaling factor (at receiver) for the first transmitter.
    Eigen::Matrix< double, 2, 3 > angularPositionScalingFactorFirstTransmitter_;

    //! Computed angular position scaling factor (at receiver) for the second transmitter.
    Eigen::Matrix< double, 2, 3 > angularPositionScalingFactorSecondTransmitter_;

    //! Computed light time correction scaling factor for first transmitter downlink.
    Eigen::Vector2d angularPositionLightTimeCorrectionScalingFirstTransmitter_;

    //! Computed light time correction scaling factor for second transmitter downlink.
    Eigen::Vector2d angularPositionLightTimeCorrectionScalingSecondTransmitter_;



    //! Relative acceleration of first transmitter w.r.t. receiver in cartesian coordinates.
    //! (retrieved from dependent variables interface)
    Eigen::Vector3d cartesianAccelerationFirstTransmitterWrtReceiver_;

    //! Relative acceleration of second transmitter w.r.t. receiver in cartesian coordinates.
    //! (retrieved from dependent variables interface)
    Eigen::Vector3d cartesianAccelerationSecondTransmitterWrtReceiver_;

    //! Partial w.r.t. link end position of relative acceleration of first transmitter w.r.t. receiver (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationFirstTransmitterWrtReceiverPosition_;

    //! Partial w.r.t. link end position of relative acceleration of first transmitter w.r.t. first transmitter (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationFirstTransmitterWrtTransmitterPosition_;

    //! Partial w.r.t. link end position of relative acceleration of first transmitter w.r.t. second transmitter (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationFirstTransmitterWrtOtherTransmitterPosition_;

    //! Partial w.r.t. link end position of relative acceleration of second transmitter w.r.t. receiver (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationSecondTransmitterWrtReceiverPosition_;

    //! Partial w.r.t. link end position of relative acceleration of second transmitter w.r.t. second transmitter (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationSecondTransmitterWrtTransmitterPosition_;

    //! Partial w.r.t. link end position of relative acceleration of second transmitter w.r.t. first transmitter (in cartesian coordinates).
    //! (retrieved from dependent variables interface)
    Eigen::Matrix3d partialAccelerationSecondTransmitterWrtOtherTransmitterPosition_;


    //! Computed right ascension of first transmitter as seen from receiver.
    double rightAscensionFirstTransmitter_;

    //! Computed right ascension of second transmitter as seen from receiver.
    double rightAscensionSecondTransmitter_;

    //! Computed declination of first transmitter as seen from receiver.
    double declinationFirstTransmitter_;

    //! Computed declination of second transmitter as seen from receiver.
    double declinationSecondTransmitter_;

    //! Computed average declination.
    double averageDeclination_;

    //! Current receiver state at central instant.
    Eigen::Vector6d receiverState_;

    //! Current first transmitter at central instant.
    Eigen::Vector6d firstTransmitterState_;

    //! Current second transmitter at central instant.
    Eigen::Vector6d secondTransmitterState_;

    //! Relative position between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d instrumentalFrameRelativePosition_;

    //! Relative velocity between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d instrumentalFrameRelativeVelocity_;

    //! Relative acceleration between the first and second transmitters in the instrumental frame of the receiver.
    Eigen::Vector2d instrumentalFrameRelativeAcceleration_;

    //! Partial of right ascension of first transmitter (as seen by receiver) w.r.t. time.
    double partialOfRightAscensionFirstTransmitterWrtTime_;

    //! Partial of declination of first transmitter (as seen by receiver) w.r.t. time.
    double partialOfDeclinationFirstTransmitterWrtTime_;

    //! Partial of right ascension of second transmitter (as seen by receiver) w.r.t. time.
    double partialOfRightAscensionSecondTransmitterWrtTime_;

    //! Partial of declination of second transmitter (as seen by receiver) w.r.t. time.
    double partialOfDeclinationSecondTransmitterWrtTime_;

    //! Second partial of right ascension of first transmitter (as seen by receiver) w.r.t. time.
    double secondPartialOfRightAscensionFirstTransmitterWrtTime_;

    //! Second partial of declination of first transmitter (as seen by receiver) w.r.t. time.
    double secondPartialOfDeclinationFirstTransmitterWrtTime_;

    //! Second partial of right ascension of second transmitter (as seen by receiver) w.r.t. time.
    double secondPartialOfRightAscensionSecondTransmitterWrtTime_;

    //! Second partial of declination of secodn transmitter (as seen by receiver) w.r.t. time.
    double secondPartialOfDeclinationSecondTransmitterWrtTime_;


    //! Partials of relative position between the two transmitters (in the instrumental frame of the receiver) w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtFirstTransmitterPosition_;

    //! Partials of relative position between the two transmitters (in the instrumental frame of the receiver) w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtSecondTransmitterPosition_;

    //! Partials of relative position between the two transmitters (in the instrumental frame of the receiver) w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativePositionWrtReceiverPosition_;


    //! Partials of relative velocity between the two transmitters (in the instrumental frame of the receiver) w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtFirstTransmitterPosition_;

    //! Partials of relative velocity between the two transmitters (in the instrumental frame of the receiver) w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtSecondTransmitterPosition_;

    //! Partials of relative velocity between the two transmitters (in the instrumental frame of the receiver) w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeVelocityWrtReceiverPosition_;


    //! Partials of relative acceleration between the two transmitters (in the instrumental frame of the receiver) w.r.t. first transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtFirstTransmitterPosition_;

    //! Partials of relative acceleration between the two transmitters (in the instrumental frame of the receiver) w.r.t. second transmitter position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtSecondTransmitterPosition_;

    //! Partials of relative acceleration between the two transmitters (in the instrumental frame of the receiver) w.r.t. receiver position.
    Eigen::Matrix< double, 2, 3 > partialsOfInstrumentalFrameRelativeAccelerationWrtReceiverPosition_;


    //! Coefficients of the cubic polynomial.
    Eigen::Vector4d cubicPolynomialCoefficients_;

    //! Coefficients of the depressed cubic polynomial.
    Eigen::Vector3d depressedCubicPolynomialCoefficients_;


    //! Partials of central instant w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtFirstTransmitterPosition_;

    //! Partials of central instant w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtSecondTransmitterPosition_;

    //! Partials of central instant w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialsOfCentralInstantWrtReceiverPosition_;

};



//! Derived class for scaling three-dimensional position partial to mutual approximation observable partial
class MutualApproximationScaling: public MutualApproximationScalingBase /*PositionPartialScaling*/
{
public:

    //! Constructor
    MutualApproximationScaling(
            const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface ) :
    MutualApproximationScalingBase( dependentVariablesInterface )
    {  }

    //! Destructor
    ~MutualApproximationScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const observation_models::LinkEnds linkEnds,
                 const Eigen::VectorXd currentObservation );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        if ( linkEndType == observation_models::transmitter )
        {
            return partialsOfCentralInstantWrtFirstTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::transmitter2 )
        {
            return partialsOfCentralInstantWrtSecondTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::receiver )
        {
            return partialsOfCentralInstantWrtReceiverPosition_;
        }
        else
        {
            throw std::runtime_error( "Error when retrieving the mutual approximation scaling factor, the link end type"
                                      "is different from transmitter, transmitter2 or receiver" );
        }
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    double getLightTimePartialScalingFactor( )
    {
//        XlightTimeCorrectionScalingFactor_ = XlightTimeCorrectionScalingFactors_.first
//                + XlightTimeCorrectionScalingFactors_.second;
//        YlightTimeCorrectionScalingFactor_ = YlightTimeCorrectionScalingFactors_.first
//                + YlightTimeCorrectionScalingFactors_.second;

        return 0.0; //1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
//                * ( XcoordinateReceiverFrame_ * XlightTimeCorrectionScalingFactor_
//                    + YcoordinateReceiverFrame_ * YlightTimeCorrectionScalingFactor_ );
    }


protected:


};


//! Derived class for scaling three-dimensional position partial to mutual approximation observable partial
class MutualApproximationWithImpactParameterScaling: public MutualApproximationScalingBase
{
public:

    //! Constructor
    MutualApproximationWithImpactParameterScaling(
            const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface ) :
    MutualApproximationScalingBase( dependentVariablesInterface )
    {  }

    //! Destructor
    ~MutualApproximationWithImpactParameterScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const observation_models::LinkEnds linkEnds,
                 const Eigen::VectorXd currentObservation );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 2, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {

        Eigen::Matrix< double, 2, 3 > scalingFactor = Eigen::Matrix< double, 2, 3 >::Zero( );

        if ( linkEndType == observation_models::transmitter )
        {
            scalingFactor.block( 0, 0, 1, 3 ) = partialsOfCentralInstantWrtFirstTransmitterPosition_;
            scalingFactor.block( 1, 0, 1, 3 ) = partialsOfImpactParameterWrtFirstTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::transmitter2 )
        {
            scalingFactor.block( 0, 0, 1, 3 ) = partialsOfCentralInstantWrtSecondTransmitterPosition_;
            scalingFactor.block( 1, 0, 1, 3 ) = partialsOfImpactParameterWrtSecondTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::receiver )
        {
            scalingFactor.block( 0, 0, 1, 3 ) = partialsOfCentralInstantWrtReceiverPosition_;
            scalingFactor.block( 1, 0, 1, 3 ) = partialsOfImpactParameterWrtReceiverPosition_;
        }
        else
        {
            throw std::runtime_error( "Error when retrieving the mutual approximation scaling factor, the link end type"
                                      "is different from transmitter, transmitter2 or receiver" );
        }

        return scalingFactor;
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    double getLightTimePartialScalingFactor( )
    {
//        XlightTimeCorrectionScalingFactor_ = XlightTimeCorrectionScalingFactors_.first
//                + XlightTimeCorrectionScalingFactors_.second;
//        YlightTimeCorrectionScalingFactor_ = YlightTimeCorrectionScalingFactors_.first
//                + YlightTimeCorrectionScalingFactors_.second;

        return 0.0; //1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
//                * ( XcoordinateReceiverFrame_ * XlightTimeCorrectionScalingFactor_
//                    + YcoordinateReceiverFrame_ * YlightTimeCorrectionScalingFactor_ );
    }


protected:

    void computePartialsOfApparentDistanceWrtLinkEndPosition( );

    void computePartialsOfImpactParameterWrtLinkEndPosition( );

private:

    //! Partials of apparent distance w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfApparentDistanceWrtFirstTransmitterPosition_;

    //! Partials of apparent distance w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfApparentDistanceWrtSecondTransmitterPosition_;

    //! Partials of apparent distance w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialsOfApparentDistanceWrtReceiverPosition_;

    //! Partials of impact parameter w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfImpactParameterWrtFirstTransmitterPosition_;

    //! Partials of impact parameter w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfImpactParameterWrtSecondTransmitterPosition_;

    //! Partials of impact parameter w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialsOfImpactParameterWrtReceiverPosition_;


};


//! Derived class for scaling three-dimensional position partial to mutual approximation observable partial
class ModifiedMutualApproximationScaling: public MutualApproximationScalingBase
{
public:

    //! Constructor
    ModifiedMutualApproximationScaling(
            const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface ) :
    MutualApproximationScalingBase( dependentVariablesInterface )
    {  }

    //! Destructor
    ~ModifiedMutualApproximationScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const observation_models::LinkEnds linkEnds,
                 const Eigen::VectorXd currentObservation );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        if ( linkEndType == observation_models::transmitter )
        {
            return partialsOfCentralInstantWrtFirstTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::transmitter2 )
        {
            return partialsOfCentralInstantWrtSecondTransmitterPosition_;
        }
        else if ( linkEndType == observation_models::receiver )
        {
            return partialsOfCentralInstantWrtReceiverPosition_;
        }
        else
        {
            throw std::runtime_error( "Error when retrieving the mutual approximation scaling factor, the link end type"
                                      "is different from transmitter, transmitter2 or receiver" );
        }
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    double getLightTimePartialScalingFactor( )
    {
//        XlightTimeCorrectionScalingFactor_ = XlightTimeCorrectionScalingFactors_.first
//                + XlightTimeCorrectionScalingFactors_.second;
//        YlightTimeCorrectionScalingFactor_ = YlightTimeCorrectionScalingFactors_.first
//                + YlightTimeCorrectionScalingFactors_.second;

        return 0.0; //1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
//                * ( XcoordinateReceiverFrame_ * XlightTimeCorrectionScalingFactor_
//                    + YcoordinateReceiverFrame_ * YlightTimeCorrectionScalingFactor_ );
    }


protected:

    void computePartialsOfModifiedObservableWrtLinkEndPosition( );

private:

    //! Partials of modified mutual approximation observable w.r.t. first transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfModifiedObservableWrtFirstTransmitterPosition_;

    //! Partials of modified mutual approximation observable w.r.t. second transmitter position.
    Eigen::Matrix< double, 1, 3 > partialsOfModifiedObservableWrtSecondTransmitterPosition_;

    //! Partials of modified mutual approximation observable w.r.t. receiver position.
    Eigen::Matrix< double, 1, 3 > partialsOfModifiedObservableWrtReceiverPosition_;

};


//! Class to compute the partial derivatives of a mutual approximation observation partial.
class MutualApproximationPartial: public ObservationPartial< 1 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > MutualApproximationPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param apparentDistanceScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    MutualApproximationPartial(
            const std::shared_ptr< MutualApproximationScalingBase > mutualApproximationScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >& lighTimeCorrectionPartials =
            std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ),
        mutualApproximationScaler_( mutualApproximationScaler ), positionPartialList_( positionPartialList ) //,
//        lighTimeCorrectionPartials_( lighTimeCorrectionPartials )
    {
        std::cout << "TEST -> CONSTRUCTOR MUTUAL APPROXIMATION PARTIAL" << "\n\n";

        // Check if the correct scaling object is provided as input.
        if ( ( std::dynamic_pointer_cast< MutualApproximationScaling >( mutualApproximationScaler_ ) == nullptr )
             && ( std::dynamic_pointer_cast< ModifiedMutualApproximationScaling >( mutualApproximationScaler_ ) == nullptr ) )
        {
            throw std::runtime_error( "Error when calculating mutual approximation partials, inconsistent scaling object." );
        }

        // Process light-time correction partials.
        if ( lighTimeCorrectionPartials.size( ) != 0 )
        {
            if ( lighTimeCorrectionPartials.size( ) != 2 )
            {
                throw std::runtime_error( "Error when making mutual approximation partials, light time corrections for "
                                          + std::to_string( lighTimeCorrectionPartials.size( ) ) + " links found, instead of 2.");
            }
            else
            {
                lighTimeCorrectionPartialsFirstTransmitter_ = lighTimeCorrectionPartials[ 0 ];
                lighTimeCorrectionPartialsSecondTransmitter_ = lighTimeCorrectionPartials[ 1 ];

                if ( lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
                     lighTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
                {
                    throw std::runtime_error( "Error when making mutual approximation partials, number of  light time "
                                              "corrections partials functions not consistent between the first transmitter leg"
                                              "and the second transmitter leg." );
                }

                // Create light time correction partial functions for link between receiver and first transmitter.
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialFirstTransmitter;
                for( unsigned int i = 0; i < lighTimeCorrectionPartialsFirstTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialFirstTransmitter = getLightTimeParameterPartialFunction(
                                parameterIdentifier, lighTimeCorrectionPartialsFirstTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialFirstTransmitter.second != 0 )
                    {
                        lighTimeCorrectionPartialsFunctionsFirstTransmitter_.push_back( lightTimeCorrectionPartialFirstTransmitter.first );
                    }
                }


                // Create light time correction partial functions for link between receiver and second transmitter.
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialSecondTransmitter;
                for( unsigned int i = 0; i < lighTimeCorrectionPartialsSecondTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialSecondTransmitter = getLightTimeParameterPartialFunction(
                                parameterIdentifier, lighTimeCorrectionPartialsSecondTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialSecondTransmitter.second != 0 )
                    {
                        lighTimeCorrectionPartialsFunctionsSecondTransmitter_.push_back( lightTimeCorrectionPartialSecondTransmitter.first );
                    }
                }

            }
        }
    }

    //! Destructor.
    ~MutualApproximationPartial( ){ }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    MutualApproximationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

    //! Function to get the number of light-time correction partial functions.
    /*!
     * Number of light-time correction partial functions.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
//        if ( lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
//             lighTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
//        {
//            throw std::runtime_error( "Error when making apparent distance partials, number of  light time "
//                                      "corrections partials functions not consistent between the first transmitter leg"
//                                      "and the second transmitter leg." );
//        }
        return lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< MutualApproximationScalingBase > mutualApproximationScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< std::function< SingleOneWayRangePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctionsFirstTransmitter_;

    std::vector< std::function< SingleOneWayRangePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctionsSecondTransmitter_;

    //! List of light-time correction partial objects.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartialsFirstTransmitter_;

    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartialsSecondTransmitter_;

    //! Pre-declare partial for current link end.
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialFirstTransmitter_;

    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialSecondTransmitter_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};



//! Class to compute the partial derivatives of a mutual approximation with impact parameter observation partial.
class MutualApproximationWithImpactParameterPartial: public ObservationPartial< 2 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > MutualApproximationWithImpactParameterPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param apparentDistanceScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    MutualApproximationWithImpactParameterPartial(
            const std::shared_ptr< MutualApproximationWithImpactParameterScaling > mutualApproximationScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >& lighTimeCorrectionPartials =
            std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) ):
        ObservationPartial< 2 >( parameterIdentifier ),
        mutualApproximationScaler_( mutualApproximationScaler ), positionPartialList_( positionPartialList ) //,
//        lighTimeCorrectionPartials_( lighTimeCorrectionPartials )
    {
        std::cout << "TEST -> CONSTRUCTOR MUTUAL APPROXIMATION WITH IMPACT PARAMETER PARTIAL" << "\n\n";

        if ( lighTimeCorrectionPartials.size( ) != 0 )
        {
            if ( lighTimeCorrectionPartials.size( ) != 2 )
            {
                throw std::runtime_error( "Error when making mutual approximation with impact parameter partials, light time corrections for "
                                          + std::to_string( lighTimeCorrectionPartials.size( ) ) + " links found, instead of 2.");
            }
            else
            {
                lighTimeCorrectionPartialsFirstTransmitter_ = lighTimeCorrectionPartials[ 0 ];
                lighTimeCorrectionPartialsSecondTransmitter_ = lighTimeCorrectionPartials[ 1 ];

                if ( lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
                     lighTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
                {
                    throw std::runtime_error( "Error when making mutual approximation with impact parameter partials, number of  light time "
                                              "corrections partials functions not consistent between the first transmitter leg"
                                              "and the second transmitter leg." );
                }

                // Create light time correction partial functions for link between receiver and first transmitter.
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialFirstTransmitter;
                for( unsigned int i = 0; i < lighTimeCorrectionPartialsFirstTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialFirstTransmitter = getLightTimeParameterPartialFunction(
                                parameterIdentifier, lighTimeCorrectionPartialsFirstTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialFirstTransmitter.second != 0 )
                    {
                        lighTimeCorrectionPartialsFunctionsFirstTransmitter_.push_back( lightTimeCorrectionPartialFirstTransmitter.first );
                    }
                }


                // Create light time correction partial functions for link between receiver and second transmitter.
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialSecondTransmitter;
                for( unsigned int i = 0; i < lighTimeCorrectionPartialsSecondTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialSecondTransmitter = getLightTimeParameterPartialFunction(
                                parameterIdentifier, lighTimeCorrectionPartialsSecondTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialSecondTransmitter.second != 0 )
                    {
                        lighTimeCorrectionPartialsFunctionsSecondTransmitter_.push_back( lightTimeCorrectionPartialSecondTransmitter.first );
                    }
                }

            }
        }
    }

    //! Destructor.
    ~MutualApproximationWithImpactParameterPartial( ){ }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    MutualApproximationWithImpactParameterPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector2d& currentObservation = Eigen::Vector2d::Constant( TUDAT_NAN ) );

    //! Function to get the number of light-time correction partial functions.
    /*!
     * Number of light-time correction partial functions.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
//        if ( lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
//             lighTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
//        {
//            throw std::runtime_error( "Error when making apparent distance partials, number of  light time "
//                                      "corrections partials functions not consistent between the first transmitter leg"
//                                      "and the second transmitter leg." );
//        }
        return lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< MutualApproximationWithImpactParameterScaling > mutualApproximationScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< std::function< SingleOneWayRangePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctionsFirstTransmitter_;

    std::vector< std::function< SingleOneWayRangePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctionsSecondTransmitter_;

    //! List of light-time correction partial objects.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartialsFirstTransmitter_;

    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartialsSecondTransmitter_;

    //! Pre-declare partial for current link end.
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialFirstTransmitter_;

    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialSecondTransmitter_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};


}

}
#endif // TUDAT_MUTUALAPPROXIMATIONPARTIAL_H
