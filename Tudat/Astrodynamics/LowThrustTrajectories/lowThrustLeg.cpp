/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"

namespace tudat
{

namespace low_thrust_trajectories
{

//! Retrieve acceleration model (thrust).
std::shared_ptr< propulsion::ThrustAcceleration > LowThrustLeg::getLowThrustAccelerationModel(
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap_[ bodyToPropagate_ ];

    // Define thrust magnitude function from the shaped trajectory.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {

        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current acceleration.
        double currentAcceleration = computeCurrentThrustAccelerationMagnitude( currentIndependentVariable, specificImpulseFunction, integratorSettings );

        // Compute current mass of the vehicle.
        double currentMass = vehicle->getBodyMass();

        // Compute and return magnitude of the low-thrust force.
        return currentAcceleration * currentMass;
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );


    // Define thrust direction function from the shaped trajectory.
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current direction of the acceleration vector.
        Eigen::Vector3d currentAccelerationDirection = computeCurrentThrustAccelerationDirection(
                    currentIndependentVariable, specificImpulseFunction, integratorSettings );

        // Return direction of the low-thrust acceleration.
        return currentAccelerationDirection;
    };

    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap_, bodyToPropagate_ );

    return lowThrustAccelerationModel;

}


//! Compute current mass of the spacecraft between two epochs.
double LowThrustLeg::computeCurrentMass(
        const double timeInitialEpoch,
        const double timeFinalEpoch,
        const double massInitialEpoch,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( massInitialEpoch );

    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = retrieveLowThrustAccelerationMap( specificImpulseFunction, integratorSettings );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Define mass propagator settings.
    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Vector1d() << massInitialEpoch ).finished(),
                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeFinalEpoch, true ) );

    // Re-initialise integrator settings.
    integratorSettings->initialTime_ = timeInitialEpoch;
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Create dynamics simulation object.
    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodyMap_, integratorSettings, massPropagatorSettings, true, false, false );

    // Propagate spacecraft mass.
    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    double currentMass = propagatedMass.rbegin()->second[ 0 ];

    return currentMass;

}

//! Compute current mass of the spacecraft.
double LowThrustLeg::computeCurrentMass( const double independentVariable,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );
    double mass = computeCurrentMass( 0.0, independentVariable, initialMass_, specificImpulseFunction, integratorSettings );
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );
    return mass;
}

//! Return mass profile.
void LowThrustLeg::getMassProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& massProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    massProfile.clear( );

    double currentMass = initialMass_;

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the mass profile of a low-thrust trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        if ( i == 0 )
        {
            currentMass = computeCurrentMass( 0.0, epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished( );
        }
        else
        {
            currentMass = computeCurrentMass( epochsVector[ i - 1 ], epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished();
        }
    }

}

//! Return thrust profile.
void LowThrustLeg::getThrustProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings)
{
    thrustProfile.clear( );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustVector = computeCurrentThrust( epochsVector[ i ], specificImpulseFunction, integratorSettings );
        thrustProfile[ epochsVector[ i ] ] = currentThrustVector;

    }
}


//! Return thrust acceleration profile.
void LowThrustLeg::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustAccelerationProfile.clear();

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustAccelerationVector = computeCurrentThrustAcceleration( epochsVector[ i ], specificImpulseFunction, integratorSettings );
        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

    }

}


void LowThrustLeg::computeSemiAnalyticalAndFullPropagation(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
                std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
        std::map< double, Eigen::VectorXd >& fullPropagationResults,
        std::map< double, Eigen::Vector6d >& semiAnalyticalResults,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory ){

    fullPropagationResults.clear();
    semiAnalyticalResults.clear();
    dependentVariablesHistory.clear();

    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ );
    integratorSettings->initialTime_ = timeOfFlight_ / 2.0;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap_, integratorSettings, propagatorSettings.first );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Declare vector of epochs at which the trajectory has to be calculated.
    std::vector< double > epochsVector;

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    // Define forward propagator settings variables.
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
    integratorSettings->initialTime_ = timeOfFlight_ / 2.0;

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap_, integratorSettings, propagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }

    getTrajectory( epochsVector, semiAnalyticalResults );

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

}


//! Define appropriate propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
std::shared_ptr< propagators::PropagatorSettings< double > > > LowThrustLeg::createLowThrustPropagatorSettings(
        std::function< double( const double ) > specificImpulseFunction,
        basic_astrodynamics::AccelerationMap perturbingAccelerationsMap,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Define propagator settings.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            propagatorSettings;

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Re-initialise spacecraft mass in body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );

    // Retrieve low-thrust trajectory nominal accelerations (thrust + central gravity accelerations).
    basic_astrodynamics::AccelerationMap lowThrustTrajectoryAccelerations = retrieveLowThrustAccelerationMap( specificImpulseFunction, integratorSettings );

    // Add perturbing accelerations given as input.
    basic_astrodynamics::AccelerationMap accelerationModelMap = perturbingAccelerationsMap;
    accelerationModelMap[ bodyToPropagate_ ][ bodyToPropagate_ ].push_back( lowThrustTrajectoryAccelerations[ bodyToPropagate_ ][ bodyToPropagate_ ][ 0 ] );
    accelerationModelMap[ bodyToPropagate_ ][ centralBody_ ].push_back( lowThrustTrajectoryAccelerations[ bodyToPropagate_ ][ centralBody_ ][ 0 ] );


    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings =
            createLowThrustTranslationalStatePropagatorSettings( accelerationModelMap, dependentVariablesToSave );

    double massHalfOfTimeOfFlight = computeCurrentMass( timeOfFlight_ / 2.0, specificImpulseFunction, integratorSettings );

    // Create settings for propagating the mass of the vehicle.
    std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
            std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate_ ] = simulation_setup::createMassRateModel(
                bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ), bodyMap_, accelerationModelMap );

    // Define backward mass propagation settings.
    massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ), terminationConditions.first );

    // Define forward mass propagation settings.
    massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ), terminationConditions.second );

    // Create list of propagation settings.
    std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
            std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

    // Backward propagator settings vector.
    propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
    propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

    // Forward propagator settings vector.
    propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
    propagatorSettingsVector.second.push_back( massPropagatorSettings.second );

    // Backward hybrid propagation settings.
    propagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
                terminationConditions.first, dependentVariablesToSave );

    // Forward hybrid propagation settings.
    propagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
                terminationConditions.second, dependentVariablesToSave );


    return propagatorSettings;
}





} // namespace low_thrust_trajectories

} // namespace tudat
