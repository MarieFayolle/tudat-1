/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEMUTUALAPPROXIMATIONPARTIALS_H
#define TUDAT_CREATEMUTUALAPPROXIMATIONPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/mutualApproximationPartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"


namespace tudat
{

namespace observation_partials
{

//! Function to generate mutual approximation partial wrt a position of a body.
/*!
 *  Function to generate mutual approximation partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param mutualApproximationLinkEnds Link ends (transmitter, transmitter2 and receiever) for which partials are to be calculated
 *  (i.e. for which mutual approximation observations are to be processed).
 *  \param bodyMap List of all bodies, for creating mutual approximation partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param mutualApproximationScaler Object scale position partials to mutual approximation partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Mutual approximation partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< MutualApproximationPartial > createMutualApproximationPartialWrtBodyPosition(
        const observation_models::LinkEnds mutualApproximationLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< MutualApproximationScaling > mutualApproximationScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ),
        const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface
                = std::shared_ptr< propagators::DependentVariablesInterface >( ) );

//! Function to generate mutual approximation partial wrt a single  parameter.
/*!
 *  Function to generate mutual approximation partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter, transmitter2 and receiever linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param mutualApproximationLinkEnds Link ends (transmitter, transmitter2 and receiever) for which mutual approximation partials are to be
 *  calculated (i.e. for which  mutual approximation observations are to be processed).
 *  \param bodyMap List of all bodies, for creating mutual approximation partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param mutualApproximationScaler Object scale position partials to mutual approximation partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Mutual approximation partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< MutualApproximationPartial > createMutualApproximationPartialWrtParameter(
        const observation_models::LinkEnds mutualApproximationLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< MutualApproximationScaling > mutualApproximationScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ),
        const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface
                = std::shared_ptr< propagators::DependentVariablesInterface >( ) )
{
    std::shared_ptr< MutualApproximationPartial > mutualApproximationPartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( mutualApproximationLinkEnds, bodyMap, parameterToEstimate );

        std::shared_ptr< MutualApproximationPartial > testMutualApproximationPartial = std::make_shared< MutualApproximationPartial >(
                    mutualApproximationScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects, dependentVariablesInterface );

        // Create mutual approximation partials if any position partials are created (i.e. if any dependency exists).
        if( positionPartials.size( ) > 0 || testMutualApproximationPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) )
        {
            mutualApproximationPartial = testMutualApproximationPartial;
        }
    }
    return mutualApproximationPartial;
}

//! Function to generate mutual approximation partials and associated scaler for single link end.
/*!
 *  Function to generate mutual approximation partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter, transmitter2 and receiver linkEndType) that are to be used.
 *  \param mutualApproximationLinkEnds Link ends (transmitter, transmitter2 and receiver) for which mutual approximation partials are to be
 *  calculated (i.e. for which mutual approximation observations are to be processed).
 *  \param bodyMap List of all bodies, for creating mutual approximation partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary mutual approximation partials of a single link end, and MutualApproximationScaling, object, used for
 *  scaling the position partial members of all MutualApproximationPartials in link end.
 */
template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
std::shared_ptr< PositionPartialScaling > >
createMutualApproximationPartials(
        const observation_models::LinkEnds mutualApproximationLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >& lightTimeCorrections =
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >( ),
        const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface
                = std::shared_ptr< propagators::DependentVariablesInterface >( ) )

{
    std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        if( lightTimeCorrections.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making mutual approximation partials, light time corrections for "
                                      + std::to_string( lightTimeCorrections.size( ) ) + " links found, instead of 2.");
        }
        lightTimeCorrectionPartialObjects.push_back(
                    observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 0 ] ) );
        lightTimeCorrectionPartialObjects.push_back(
                    observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 1 ] ) );
    }

    // Create scaling object, to be used for all mutual approximation partials in current link end.
    std::shared_ptr< MutualApproximationScaling > mutualApproximationScaling = std::make_shared< MutualApproximationScaling >( dependentVariablesInterface );

    SingleLinkObservationPartialList mutualApproximationPartials;


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
            throw std::runtime_error( "Error when making mutual approximation partials, could not identify parameter" );
        }

        // Create position mutual approximation partial for current body
        std::shared_ptr< MutualApproximationPartial > currentMutualApproximationPartial = createMutualApproximationPartialWrtBodyPosition(
                    mutualApproximationLinkEnds, bodyMap, acceleratedBody, mutualApproximationScaling,
                    lightTimeCorrectionPartialObjects, dependentVariablesInterface );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current mutual approximation and current body)
        if( currentMutualApproximationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            mutualApproximationPartials[ currentPair ] = currentMutualApproximationPartial;
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
        // Create position mutual approximation partial for current parameter
        std::shared_ptr< MutualApproximationPartial > currentMutualApproximationPartial = createMutualApproximationPartialWrtParameter(
                    mutualApproximationLinkEnds, bodyMap, parameterIterator->second, mutualApproximationScaling,
                    lightTimeCorrectionPartialObjects, dependentVariablesInterface );

        if( currentMutualApproximationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            mutualApproximationPartials[ currentPair ] = currentMutualApproximationPartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate = parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator = vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        // Create position mutual approximation partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentMutualApproximationPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentMutualApproximationPartial = createMutualApproximationPartialWrtParameter(
                        mutualApproximationLinkEnds, bodyMap, parameterIterator->second, mutualApproximationScaling,
                        lightTimeCorrectionPartialObjects, dependentVariablesInterface );
        }
        else
        {
            currentMutualApproximationPartial = createObservationPartialWrtLinkProperty< 1 >(
                        mutualApproximationLinkEnds, observation_models::mutual_approximation, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentMutualApproximationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            mutualApproximationPartials[ currentPair ] = currentMutualApproximationPartial;
        }

    }
    return std::make_pair( mutualApproximationPartials, mutualApproximationScaling );
}

//! Function to generate mutual approximation partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate mutual approximation partials for all parameters that are to be estimated, for all sets of link ends.
 *  The mutual approximation partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter, transmitter2 and receiever linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which mutual approximation distance partials are to be calculated (i.e. for which one-way
 *  mutual approximation observations are  to be processed).
 *  \param bodyMap List of all bodies, for creating mutual approximation partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary mutual approximation partials of a single link
 *  end, and MutualApproximationScaling, object, used for scaling the position partial members of all MutualApproximationPartials in
 *  link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
std::shared_ptr< PositionPartialScaling > > >
createMutualApproximationPartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ),
        const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface
                = std::shared_ptr< propagators::DependentVariablesInterface >( ) )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            std::shared_ptr< PositionPartialScaling > > > mutualApproximationPartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::receiver ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter2 ) == 0 ) )
        {
            throw std::runtime_error( "Error when making mutual approximation partials, did not find both transmitter, transmitter2 and receiver in link ends" );

        }

        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) != 2 )
            {
                std::cerr << "Error when making mutual approximation partials, light time corrections for " <<
                           lightTimeCorrections.at( linkEnds.at( i ) ).size( ) << " links found, instead of 2." << std::endl;
            }
            currentLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) );
        }
        else
        {
            currentLightTimeCorrections.clear( );
        }

        // Create mutual approximation partials for current link ends
        mutualApproximationPartials[ linkEnds[ i ] ] = createMutualApproximationPartials(
                    linkEnds[ i ], bodyMap, parametersToEstimate, currentLightTimeCorrections, dependentVariablesInterface );
    }
    return mutualApproximationPartials;
}

}

}

#endif // TUDAT_CREATEMUTUALAPPROXIMATIONPARTIALS_H
