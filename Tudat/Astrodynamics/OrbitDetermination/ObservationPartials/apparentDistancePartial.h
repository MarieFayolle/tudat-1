/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_APPARENTDISTANCEPARTIAL_H
#define TUDAT_APPARENTDISTANCEPARTIAL_H

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/angularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the right ascension and declination of the observed body as seen from the observer.
/*!
 * Function to compute the right ascension and declination of the observed body as seen from the observer.
 * \param relativePosition Vector from observer to observed object
 * \return pair object containing the right ascension and declination of the observed body as seen from the observer.
 */
std::pair< double, double > computeRightAscensionAndDeclination( Eigen::Vector3d relativePosition );


//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
/*!
 * Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > computePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver );

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
/*!
 * Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return derivative of (direct geometric) declination w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > computePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

//! Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
//! observed object.
/*!
 * Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
 * observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return Derivative of (direct geometric) right ascension and declination w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > computePartialOfApparentDistanceWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

//! Derived class for scaling three-dimensional position partial to apparent distance observable partial
class ApparentDistanceScaling: public PositionPartialScaling
{
public:

    //! Constructor
    ApparentDistanceScaling( )
    {
        referenceAngularPositionScaler_ = std::make_shared< AngularPositionScaling >( );
    }

    //! Destructor
    ~ApparentDistanceScaling( ){ }

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
            XscalingFactor_ = - XscalingFactors_.first;
            YscalingFactor_ = - YscalingFactors_.first;
//            return angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) * std::cos( averageDeclination_ )
//                    + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
//                    * std::sin( averageDeclination_ )
//                    * angularPositionLightTimeCorrectionScalingFirstTransmitter_.block( 1, 0, 1, 3 );
        }
        else if ( linkEndType == observation_models::transmitter2 )
        {
            XscalingFactor_ = - XscalingFactors_.second;
            YscalingFactor_ = - YscalingFactors_.second;
//            return - angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 ) * std::cos( averageDeclination_ )
//                    + ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
//                    * std::sin( averageDeclination_ )
//                    * angularPositionLightTimeCorrectionScalingSecondTransmitter_.block( 1, 0, 1, 3 );
        }
        else if ( linkEndType == observation_models::receiver )
        {
            XscalingFactor_ = XscalingFactors_.first + XscalingFactors_.second;
            YscalingFactor_ = YscalingFactors_.first + YscalingFactors_.second;
//            return ( angularPositionScalingFactorSecondTransmitter_.block( 0, 0, 1, 3 )
//                     - angularPositionScalingFactorFirstTransmitter_.block( 0, 0, 1, 3 ) ) * std::cos( averageDeclination_ )
//                    - ( rightAscensionSecondTransmitter_ - rightAscensionFirstTransmitter_ )
//                    * std::sin( averageDeclination_ )
//                    * ( angularPositionLightTimeCorrectionScalingFirstTransmitter_.block( 1, 0, 1, 3 )
//                        + angularPositionLightTimeCorrectionScalingSecondTransmitter_.block( 1, 0, 1, 3 ) );
        }
//        else
//        {
//            throw std::runtime_error( "Error when retrieving the apparent distance scaling factor, the link end type"
//                                      "is different from transmitter, transmitter2 or receiver" );
//        }

        return 1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
                * ( XcoordinateReceiverFrame_ * XscalingFactor_ + YcoordinateReceiverFrame_ * YscalingFactor_ );
//        if ( fixedLinkEnd != observation_models::receiver )
//        {
//            throw std::runtime_error( "Error when updating the apparent distance scaling object, fixed link end time different from receiver." );
//        }

//        return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    double getLightTimePartialScalingFactor( )
    {
        XlightTimeCorrectionScalingFactor_ = XlightTimeCorrectionScalingFactors_.first
                + XlightTimeCorrectionScalingFactors_.second;
        YlightTimeCorrectionScalingFactor_ = YlightTimeCorrectionScalingFactors_.first
                + YlightTimeCorrectionScalingFactors_.second;

        return 1.0 / std::sqrt( XcoordinateReceiverFrame_ * XcoordinateReceiverFrame_ + YcoordinateReceiverFrame_ * YcoordinateReceiverFrame_ )
                * ( XcoordinateReceiverFrame_ * XlightTimeCorrectionScalingFactor_
                    + YcoordinateReceiverFrame_ * YlightTimeCorrectionScalingFactor_ );
    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }


private:

    //! Computed angular position scaling factor (at receiver) for the first transmitter.
    Eigen::Matrix< double, 2, 3 > angularPositionScalingFactorFirstTransmitter_;

    //! Computed angular position scaling factor (at receiver) for the second transmitter.
    Eigen::Matrix< double, 2, 3 > angularPositionScalingFactorSecondTransmitter_;

//    //! Predeclared common scaling factor
//    Eigen::Matrix< double, 1, 3 > scalingFactor_;

    //! Computed scaling factor for first transmitter (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorFirstTransmitter_;

    //! Computed scaling factor for second transmitter (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorSecondTransmitter_;

//    //! Computed scaling factor (at receiver)
//    Eigen::Matrix< double, 1, 3 > referenceScalingFactor_;

//    //! Computed light time correction scaling factor
//    Eigen::Vector1d referenceLightTimeCorrectionScaling_;

    //! Computed light time correction scaling factor for first transmitter downlink.
    Eigen::Vector2d angularPositionLightTimeCorrectionScalingFirstTransmitter_;

    //! Computed light time correction scaling factor for second transmitter downlink.
    Eigen::Vector2d angularPositionLightTimeCorrectionScalingSecondTransmitter_;

    //! Pair containing the partials of the right ascension contribution to the coordinate X (in instrumental frame).
    //! The first element of the pair refers to the first transmitter - receiver segment, while the second element refers
    //! to the second transmitter - receiver segment.
    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > XscalingFactors_;

//    //! Pair containing the partials of the declination contribution to the coordinate X (in instrumental frame).
//    //! The first element of the pair refers to the first transmitter - receiver segment, while the second element refers
//    //! to the second transmitter - receiver segment.
//    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > scalingFactorsXdeclinationContribution_;

    //! Pair containing the partials of the coordinate Y (in instrumental frame).
    //! The first element of the pair refers to the first transmitter - receiver segment, while the second element refers
    //! to the second transmitter - receiver segment.
    std::pair< Eigen::Matrix< double, 1, 3 >, Eigen::Matrix< double, 1, 3 > > YscalingFactors_;


    std::pair< double, double > XlightTimeCorrectionScalingFactors_;
    std::pair< double, double > YlightTimeCorrectionScalingFactors_;

    std::pair< double, double > scalingFactorXCoefficient_;
    double scalingFactorYCoefficient_;

    Eigen::Matrix< double, 1, 3 > XscalingFactor_;
    Eigen::Matrix< double, 1, 3 > YscalingFactor_;

    double XlightTimeCorrectionScalingFactor_;
    double YlightTimeCorrectionScalingFactor_;



    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

    //!
    std::shared_ptr< AngularPositionScaling > referenceAngularPositionScaler_;

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

    //! Computed X coordinate in instrumental frame.
    double XcoordinateReceiverFrame_;

    //! Computed Y coordinate in instrumental frame.
    double YcoordinateReceiverFrame_;



};

//! Class to compute the partial derivatives of an apparent distance observation partial.
class ApparentDistancePartial: public ObservationPartial< 1 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > ApparentDistancePartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param apparentDistanceScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    ApparentDistancePartial(
            const std::shared_ptr< ApparentDistanceScaling > apparentDistanceScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >& lighTimeCorrectionPartials =
            std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ),
        apparentDistanceScaler_( apparentDistanceScaler ), positionPartialList_( positionPartialList ) //,
//        lighTimeCorrectionPartials_( lighTimeCorrectionPartials )
    {
        if ( lighTimeCorrectionPartials.size( ) != 0 )
        {
            if ( lighTimeCorrectionPartials.size( ) != 2 )
            {
                throw std::runtime_error( "Error when making apparent distance partials, light time corrections for "
                                          + std::to_string( lighTimeCorrectionPartials.size( ) ) + " links found, instead of 2.");
            }
            else
            {
                lighTimeCorrectionPartialsFirstTransmitter_ = lighTimeCorrectionPartials[ 0 ];
                lighTimeCorrectionPartialsSecondTransmitter_ = lighTimeCorrectionPartials[ 1 ];

                if ( lighTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
                     lighTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
                {
                    throw std::runtime_error( "Error when making apparent distance partials, number of  light time "
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

//            std::pair< std::function< SingleOneWayRangePartialReturnType(
//                        const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >,
//                    bool > lightTimeCorrectionPartial;

//            // Create light time correction partial functions
//            for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
//            {
//                lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
//                            parameterIdentifier, lighTimeCorrectionPartials.at( i ) );
//                if( lightTimeCorrectionPartial.second != 0 )
//                {
//                    lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
//                }
//            }
        }
    }

    //! Destructor.
    ~ApparentDistancePartial( ){ }

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
    ApparentDistancePartialReturnType calculatePartial(
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
    std::shared_ptr< ApparentDistanceScaling > apparentDistanceScaler_;

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
#endif // TUDAT_APPARENTDISTANCEPARTIAL_H
