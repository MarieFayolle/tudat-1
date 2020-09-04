/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EstimationSetup/createMutualApproximationPartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"

namespace tudat
{


namespace observation_partials
{

//! Function to generate mutual approximation partial wrt a position of a body.
std::shared_ptr< MutualApproximationPartial > createMutualApproximationPartialWrtBodyPosition(
        const observation_models::LinkEnds mutualApproximationLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< MutualApproximationScaling > mutualApproximationScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects,
        const std::shared_ptr< propagators::DependentVariablesInterface > dependentVariablesInterface )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( mutualApproximationLinkEnds, bodyMap, bodyToEstimate );

    // Create mutual approximation if any position partials are created (i.e. if any dependency exists).
    std::shared_ptr< MutualApproximationPartial > mutualApproximationPartial;
    if( positionPartials.size( ) > 0 )
    {
        mutualApproximationPartial = std::make_shared< MutualApproximationPartial >(
                    mutualApproximationScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects, dependentVariablesInterface );
    }

    return mutualApproximationPartial;
}

}

}

