/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

//template class VariationalEquationsSolver< double, double >;
//template class VariationalEquationsSolver< long double, double >;
//template class VariationalEquationsSolver< double, Time >;
//template class VariationalEquationsSolver< long double, Time >;

//template class SingleArcVariationalEquationsSolver< double, double >;
//template class SingleArcVariationalEquationsSolver< long double, double >;
//template class SingleArcVariationalEquationsSolver< double, Time >;
//template class SingleArcVariationalEquationsSolver< long double, Time >;

//template class MultiArcVariationalEquationsSolver< double, double >;
//template class MultiArcVariationalEquationsSolver< long double, double >;
//template class MultiArcVariationalEquationsSolver< double, Time >;
//template class MultiArcVariationalEquationsSolver< long double, Time >;

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
void createStateTransitionAndSensitivityMatrixInterpolator(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const bool clearRawSolution )
{
    // Create interpolator for state transition matrix.
    stateTransitionMatrixInterpolator=
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( variationalEquationsSolution[ 0 ] ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( variationalEquationsSolution[ 0 ] ), 4 );
    if( clearRawSolution )
    {
        variationalEquationsSolution[ 0 ].clear( );
    }

//    std::cout<<"State trans. size "<<variationalEquationsSolution[ 0 ].size( )<<std::endl;
//    std::cout<<"State trans. matrix "<<variationalEquationsSolution[ 0 ].begin( )->second<<std::endl;

    // Create interpolator for sensitivity matrix.
    sensitivityMatrixInterpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( variationalEquationsSolution[ 1 ] ),
                utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( variationalEquationsSolution[ 1 ] ), 4 );

    //std::cout<<"State trans "<<stateTransitionMatrixInterpolator->interpolate( 20000.0 )<<std::endl;

    if( clearRawSolution )
    {
        variationalEquationsSolution[ 1 ].clear( );
    }

}



void getDependentVariablesHistoryForInterface(
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistoryInterface,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesSaveSettings,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesInterfaceSettings )
{
    dependentVariablesHistoryInterface.clear( );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > singleDependentVariablesSaveSettings = dependentVariablesSaveSettings->dependentVariables_;
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > singleDependentVariablesInterfaceSettings = dependentVariablesInterfaceSettings->dependentVariables_;

//        std::map< PropagationDependentVariables, int > dependentVariableInterfaceIndices;
    std::map< std::string, int > dependentVariablesSaveIndices;
    std::vector< std::pair< int, int > > dependentVariableInterfaceIndices;

    int totalDependentVariablesSize = 0;
    for ( unsigned int i = 0 ; i < singleDependentVariablesSaveSettings.size( ) ; i++ )
    {
        dependentVariablesSaveIndices[ getDependentVariableId( singleDependentVariablesSaveSettings[ i ] ) ]
                = totalDependentVariablesSize;
        totalDependentVariablesSize += getDependentVariableSaveSize( singleDependentVariablesSaveSettings[ i ] );
    }

    int totalDependentVariablesInterfaceSize = 0;
    for ( unsigned int i = 0 ; i < singleDependentVariablesInterfaceSettings.size( ) ; i++ )
    {
        std::string currentInterfaceDependentVariableId =
                getDependentVariableId( singleDependentVariablesInterfaceSettings[ i ] );
        if ( dependentVariablesSaveIndices.count( currentInterfaceDependentVariableId ) != 1 )
        {
            throw std::runtime_error( "Error when creating a dependent variables interface, required dependent variables settings"
                                      "incompatible with dependent variables to save." );
        }

        std::map< std::string, int >::iterator itr = dependentVariablesSaveIndices.find( currentInterfaceDependentVariableId );
//            dependentVariableInterfaceIndices[ currentInterfaceDependentVariableType ] = itr->second;
        dependentVariableInterfaceIndices.push_back(
                    std::make_pair( itr->second,
                                    getDependentVariableSaveSize( singleDependentVariablesInterfaceSettings[ i ] ) ) );
        totalDependentVariablesInterfaceSize += getDependentVariableSaveSize( singleDependentVariablesInterfaceSettings[ i ] );
    }

//        for ( std::map< PropagationDependentVariables, int >::iterator itr =
//              dependentVariableInterfaceIndices.begin( ) ; itr != dependentVariableInterfaceIndices.end( ) ; itr++ )
//        {

//        }
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ;
          itr != dependentVariablesHistory.end( ) ; itr++ )
    {
        Eigen::VectorXd currentDependentVariableValues = Eigen::VectorXd::Zero( totalDependentVariablesInterfaceSize );

        int currentIndex = 0;
        for ( unsigned int i = 0 ; i < dependentVariableInterfaceIndices.size( ) ; i++ )
        {
            currentDependentVariableValues.segment( currentIndex, dependentVariableInterfaceIndices[ i ].second )
                    = itr->second.segment( dependentVariableInterfaceIndices[ i ].first,
                                           dependentVariableInterfaceIndices[ i ].second );
            currentIndex += dependentVariableInterfaceIndices[ i ].second;
        }

        dependentVariablesHistoryInterface[ itr->first ] = currentDependentVariableValues;
    }

}


//! Function to create interpolator for dependent variable from numerical results.
void createDependentVariablesInterpolator(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >&
        dependentVariablesInterpolator,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory,
        const bool clearDependentVariablesHistory )
{
    // Create interpolator for dependent variable.
    dependentVariablesInterpolator=
            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory ),
                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory ), 4 );

    if( clearDependentVariablesHistory )
    {
        dependentVariablesHistory.clear( );
    }
}


template class VariationalEquationsSolver< double, double >;
template class SingleArcVariationalEquationsSolver< double, double >;
template class MultiArcVariationalEquationsSolver< double, double >;
template class HybridArcVariationalEquationsSolver< double, double >;

#if( BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class VariationalEquationsSolver< long double, double >;
template class VariationalEquationsSolver< double, Time >;
template class VariationalEquationsSolver< long double, Time >;

template class SingleArcVariationalEquationsSolver< long double, double >;
template class SingleArcVariationalEquationsSolver< double, Time >;
template class SingleArcVariationalEquationsSolver< long double, Time >;

template class MultiArcVariationalEquationsSolver< long double, double >;
template class MultiArcVariationalEquationsSolver< double, Time >;
template class MultiArcVariationalEquationsSolver< long double, Time >;

template class HybridArcVariationalEquationsSolver< long double, double >;
template class HybridArcVariationalEquationsSolver< double, Time >;
template class HybridArcVariationalEquationsSolver< long double, Time >;
#endif



}

}
