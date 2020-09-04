/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEAPPARENTDISTANCEPARTIALS_H
#define TUDAT_CREATEAPPARENTDISTANCEPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/apparentDistancePartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"


namespace tudat
{

namespace observation_partials
{

//! Function to generate apparent distance partial wrt a position of a body.
/*!
 *  Function to generate apparent distance partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param apparentDistanceLinkEnds Link ends (transmitter, transmitter2 and receiever) for which partials are to be calculated
 *  (i.e. for which apparent distance observations are to be processed).
 *  \param bodyMap List of all bodies, for creating apparent distance partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param apparentDistanceScaler Object scale position partials to apparent distance partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Apparent distance partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< ApparentDistancePartial > createApparentDistancePartialWrtBodyPosition(
        const observation_models::LinkEnds apparentDistanceLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< ApparentDistanceScaling > apparentDistanceScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) );

//! Function to generate apparent distance partial wrt a single  parameter.
/*!
 *  Function to generate apparent distance partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter, transmitter2 and receiever linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param apparentDistanceLinkEnds Link ends (transmitter, transmitter2 and receiever) for which apparent distance partials are to be
 *  calculated (i.e. for which  apparent distance observations are to be processed).
 *  \param bodyMap List of all bodies, for creating apparent distance partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param apparentDistanceScaler Object scale position partials to apparent distance partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Apparent distance partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< ApparentDistancePartial > createApparentDistancePartialWrtParameter(
        const observation_models::LinkEnds apparentDistanceLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< ApparentDistanceScaling > apparentDistanceScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) )
{
    std::shared_ptr< ApparentDistancePartial > apparentDistancePartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( apparentDistanceLinkEnds, bodyMap, parameterToEstimate );

        std::shared_ptr< ApparentDistancePartial > testApparentDistancePartial = std::make_shared< ApparentDistancePartial >(
                    apparentDistanceScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );

        // Create apparent distance partials if any position partials are created (i.e. if any dependency exists).
        if( positionPartials.size( ) > 0 || testApparentDistancePartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) )
        {
            apparentDistancePartial = testApparentDistancePartial;
        }
    }
    return apparentDistancePartial;
}

//! Function to generate apparent distance partials and associated scaler for single link end.
/*!
 *  Function to generate apparent distance partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter, transmitter2 and receiever linkEndType) that are to be used.
 *  \param apparentDistanceLinkEnds Link ends (transmitter, transmitter2 and receiever) for which apparent distance partials are to be
 *  calculated (i.e. for which apparent distance observations are to be processed).
 *  \param bodyMap List of all bodies, for creating apparent distance partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary apparent distance partials of a single link end, and ApparentDistanceScaling, object, used for
 *  scaling the position partial members of all ApparentDistancePartials in link end.
 */
template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
std::shared_ptr< PositionPartialScaling > >
createApparentDistancePartials(
        const observation_models::LinkEnds apparentDistanceLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >& lightTimeCorrections =
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >( ) )

{
    std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        if( lightTimeCorrections.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making apparent distance partials, light time corrections for "
                                      + std::to_string( lightTimeCorrections.size( ) ) + " links found, instead of 2.");
        }
        lightTimeCorrectionPartialObjects.push_back(
                    observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 0 ] ) );
        lightTimeCorrectionPartialObjects.push_back(
                    observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 1 ] ) );
    }

    // Create scaling object, to be used for all apparent distance partials in current link end.
    std::shared_ptr< ApparentDistanceScaling > apparentDistanceScaling = std::make_shared< ApparentDistanceScaling >( );

    SingleLinkObservationPartialList apparentDistancePartials;


    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {

        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else
        {
            throw std::runtime_error( "Error when making apparent distance partials, could not identify parameter" );
        }

        // Create position apparent distance partial for current body
        std::shared_ptr< ApparentDistancePartial > currentApparentDistancePartial = createApparentDistancePartialWrtBodyPosition(
                    apparentDistanceLinkEnds, bodyMap, acceleratedBody, apparentDistanceScaling ,
                    lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current apparent distance and current body)
        if( currentApparentDistancePartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            apparentDistancePartials[ currentPair ] = currentApparentDistancePartial;
        }

        // Increment current index by size of body initial state (6).
        currentIndex += 6;
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator = doubleParametersToEstimate.begin( );
         parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position apparent distance partial for current parameter
        std::shared_ptr< ApparentDistancePartial > currentApparentDistancePartial = createApparentDistancePartialWrtParameter(
                    apparentDistanceLinkEnds, bodyMap, parameterIterator->second, apparentDistanceScaling,
                    lightTimeCorrectionPartialObjects );

        if( currentApparentDistancePartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            apparentDistancePartials[ currentPair ] = currentApparentDistancePartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate = parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator = vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        // Create position apparent distance partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentApparentDistancePartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentApparentDistancePartial = createApparentDistancePartialWrtParameter(
                        apparentDistanceLinkEnds, bodyMap, parameterIterator->second, apparentDistanceScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentApparentDistancePartial = createObservationPartialWrtLinkProperty< 1 >(
                        apparentDistanceLinkEnds, observation_models::apparent_distance, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentApparentDistancePartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            apparentDistancePartials[ currentPair ] = currentApparentDistancePartial;
        }

    }
    return std::make_pair( apparentDistancePartials, apparentDistanceScaling );
}

//! Function to generate apparent distance partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate apparent distance partials for all parameters that are to be estimated, for all sets of link ends.
 *  The apparent distance partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter, transmitter2 and receiever linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which apparent distance partials are to be calculated (i.e. for which one-way
 *  apparent distance observations are  to be processed).
 *  \param bodyMap List of all bodies, for creating apparent distance partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary apparent distance partials of a single link
 *  end, and ApparentDistanceScaling, object, used for scaling the position partial members of all ApparentDistancePartials in
 *  link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
std::shared_ptr< PositionPartialScaling > > >
createApparentDistancePartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            std::shared_ptr< PositionPartialScaling > > > apparentDistancePartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::receiver ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter2 ) == 0 ) )
        {
            throw std::runtime_error( "Error when making apparent distance partials, did not find both transmitter, transmitter2 and receiver in link ends" );

        }

        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) != 2 )
            {
                std::cerr << "Error when making apparent distance partials, light time corrections for " <<
                           lightTimeCorrections.at( linkEnds.at( i ) ).size( ) << " links found, instead of 2." << std::endl;
            }
            currentLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) );
        }
        else
        {
            currentLightTimeCorrections.clear( );
        }

        // Create apparent distance partials for current link ends
        apparentDistancePartials[ linkEnds[ i ] ] = createApparentDistancePartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate, currentLightTimeCorrections );
    }
    return apparentDistancePartials;
}

}

}

#endif // TUDAT_CREATEAPPARENTDISTANCEPARTIALS_H
