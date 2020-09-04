/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Propagators/dependentVariablesInterface.h"

namespace tudat
{

namespace propagators
{

//! Function to reset the dependent variable interpolators
void SingleArcDependentVariablesInterface::updateDependentVariablesInterpolators(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >
        dependentVariablesInterpolator/*,
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
        sensitivityMatrixInterpolator,
        const std::vector< std::pair< int, int > >& statePartialAdditionIndices*/ )
{
    dependentVariablesInterpolator_ = dependentVariablesInterpolator;
//    statePartialAdditionIndices_ = statePartialAdditionIndices;
}

//! Function to get the dependent variable at a given time.
Eigen::VectorXd SingleArcDependentVariablesInterface::getDependentVariables(
        const double evaluationTime )
{
    dependentVariables_.setZero( );

    // Set dependent variable.
    dependentVariables_ = dependentVariablesInterpolator_->interpolate( evaluationTime );
//    combinedStateTransitionMatrix_.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
//            stateTransitionMatrixInterpolator_->interpolate( evaluationTime );

//    if( sensitivityMatrixSize_ > 0 )
//    {
//        combinedStateTransitionMatrix_.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
//                sensitivityMatrixInterpolator_->interpolate( evaluationTime );
//    }

//    for( unsigned int i = 0; i < statePartialAdditionIndices_.size( ); i++ )
//    {
//        combinedStateTransitionMatrix_.block( statePartialAdditionIndices_.at( i ).first, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ ) +=
//                combinedStateTransitionMatrix_.block( statePartialAdditionIndices_.at( i ).second, 0, 6, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );
//    }


    return dependentVariables_;
}

//! Constructor
MultiArcDependentVariablesInterface::MultiArcDependentVariablesInterface(
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
        dependentVariablesInterpolators,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesSettings,
//            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
//            sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes,
        const std::vector< double >& arcEndTimes/*,
        const int numberOfInitialDynamicalParameters*//*,
        const int numberOfParameters,
        const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices*/ ):
    DependentVariablesInterface( dependentVariablesSettings ),
    dependentVariablesInterpolators_( dependentVariablesInterpolators ),
    arcStartTimes_( arcStartTimes ),
    arcEndTimes_( arcEndTimes )/*,
    statePartialAdditionIndices_( statePartialAdditionIndices )*/
{
    if( arcStartTimes_.size( ) != arcEndTimes_.size( ) )
    {
        throw std::runtime_error( "Error when making MultiArcDependentVariablesInterface, incompatible time lists" );
    }
    numberOfStateArcs_ = arcStartTimes_.size( );

//    sensitivityMatrixSize_ = numberOfParameters - numberOfStateArcs_ * stateTransitionMatrixSize_;

//    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
//            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
//    {
//        throw std::runtime_error(
//                    "Error when making multi arc state transition and sensitivity interface, vector sizes are inconsistent" );
//    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
    lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to reset the dependent variables interpolators
void MultiArcDependentVariablesInterface::updateDependentVariablesInterpolators(
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
        dependentVariablesInterpolators,
//            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >
//            sensitivityMatrixInterpolators,
        const std::vector< double >& arcStartTimes,
        const std::vector< double >& arcEndTimes/*,
        const std::vector< std::vector< std::pair< int, int > > >& statePartialAdditionIndices*/  )
{
    dependentVariablesInterpolators_ = dependentVariablesInterpolators;
    arcStartTimes_ =  arcStartTimes;
    arcEndTimes_ = arcEndTimes;
//    statePartialAdditionIndices_ = statePartialAdditionIndices;

//    if( stateTransitionMatrixInterpolators_.size( ) != sensitivityMatrixInterpolators_.size( ) ||
//            stateTransitionMatrixInterpolators_.size( ) != static_cast< unsigned int >( numberOfStateArcs_ ) )
//    {
//        throw std::runtime_error(
//                    "Error when resetting multi arc state transition and sensitivity interface, vector sizes are inconsistent." );
//    }

    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

    lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                arcSplitTimes );
}

//! Function to get the single-arc dependent variable at a given time.
Eigen::VectorXd MultiArcDependentVariablesInterface::getDependentVariables(
        const double evaluationTime/*,
        const bool addCentralBodySensitivity*/ )
{
    Eigen::VectorXd dependentVariablesVector = Eigen::VectorXd::Zero( dependentVariablesSize_ );
//                stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

    int currentArc = getCurrentArc( evaluationTime ).first;

    // Set dependent variables vector.
    if( currentArc >= 0 )
    {
        dependentVariablesVector.segment( 0, dependentVariablesSize_) = dependentVariablesInterpolators_.at( currentArc )->interpolate( evaluationTime );
//        , stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
//                stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
//        combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
//                sensitivityMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );

//        for( unsigned int i = 0; i < statePartialAdditionIndices_.at( currentArc ).size( ); i++ )
//        {
//            int indicesToAdd = addCentralBodySensitivity ?
//                        ( stateTransitionMatrixSize_ + sensitivityMatrixSize_ ) : stateTransitionMatrixSize_;
//            combinedStateTransitionMatrix.block(
//                        statePartialAdditionIndices_.at( currentArc ).at( i ).first, 0, 6, indicesToAdd ) +=
//                    combinedStateTransitionMatrix.block(
//                        statePartialAdditionIndices_.at( currentArc ).at( i ).second, 0, 6, indicesToAdd );
//        }
    }
    return dependentVariablesVector;
}

////! Function to get the single-arc dependent variable at a given time.
//Eigen::VectorXd MultiArcDependentVariablesInterface::getDependentVariable(
//        const double evaluationTime )
//{
//    return getDependentVariable( evaluationTime, true );
//}

////! Function to get the dependent variable for each arc at a given time.
//Eigen::VectorXd MultiArcDependentVariableInterface::getFullDependentVariable(
//        const double evaluationTime,
//        const bool addCentralBodySensitivity )
//{
//    Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix(
//                evaluationTime, addCentralBodySensitivity );
//    Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
//                stateTransitionMatrixSize_, numberOfStateArcs_ * stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

//    int currentArc = getCurrentArc( evaluationTime ).first;

//    // Set Phi and S matrices of current arc.
//    if( currentArc >= 0 )
//    {
//        fullCombinedStateTransitionMatrix.block(
//                    0, currentArc * stateTransitionMatrixSize_, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ) =
//                stateTransitionMatrixInterpolators_.at( currentArc )->interpolate( evaluationTime );
//        fullCombinedStateTransitionMatrix.block(
//                    0, numberOfStateArcs_ * stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ ) =
//                combinedStateTransitionMatrix.block( 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_, sensitivityMatrixSize_ );
//    }
//    return fullCombinedStateTransitionMatrix;
//}

////! Function to get the dependent variable for each arc at a given time.
//Eigen::VectorXd MultiArcDependentVariableInterface::getFullDependentVariable(
//        const double evaluationTime )
//{
//    return getFullDependentVariable( evaluationTime, true );
//}

//! Function to retrieve the current arc for a given time
std::pair< int, double > MultiArcDependentVariablesInterface::getCurrentArc( const double evaluationTime )
{

    int currentArc =  lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
    if( evaluationTime <= arcEndTimes_.at( currentArc ) && evaluationTime >= arcStartTimes_.at( currentArc ) )
    {
        return std::make_pair( currentArc, arcStartTimes_.at( currentArc ) );
    }
    else
    {
        return std::make_pair( -1, TUDAT_NAN );
    }
}

//! Function to get the dependent variable at a given time.
Eigen::VectorXd HybridArcDependentVariablesInterface::getDependentVariables(
        const double evaluationTime )
{
    Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( hybridArcDependentVariablesSize_ );
//                stateTransitionMatrixSize_, stateTransitionMatrixSize_ + sensitivityMatrixSize_ );

    // Get single-arc dependent variables.
    Eigen::VectorXd singleArcDependentVariables = singleArcInterface_->getDependentVariables(
                evaluationTime );

    // Get multi-arc dependent variables.
    Eigen::VectorXd multiArcDependentVariables = multiArcInterface_->getDependentVariables(
                evaluationTime );
    std::pair< int, double >  currentArc = multiArcInterface_->getCurrentArc( evaluationTime );

    if ( ( /*singleArcDependentVariables.size( ) +*/ multiArcDependentVariables.size( ) ) != hybridArcDependentVariablesSize_ )
    {
        throw std::runtime_error( "Error when getting dependent variables from hybrid arc interfaces, size inconsistent "
                                  "with size of single- and multi-arc dependent variables." );
    }

    // Set single-arc block
    dependentVariables.segment( 0, singleArcDependentVariablesSize_ ) = //, singleArcStateSize_, singleArcStateSize_ ) =
            singleArcDependentVariables; //( 0, 0, singleArcStateSize_, singleArcStateSize_ );

//    // Set single-arc sensitivity block
//    dependentVariables.block( 0, multiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ ) =
//            singleArcDependentVariables.block( 0, singleArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );

    if( !( currentArc.first < 0 ) )
    {
        int currentMultiArcDependentVariableIndex = singleArcDependentVariablesSize_;

        // Set multi-arc block
        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices_.begin( ) ;
              itr != multiArcDependentVariablesIdsAndIndices_.end( ) ; itr++ )
        {
            int sizeCurrentDependentVariable = 0;
            for ( unsigned int i = 0 ;
                  i < multiArcInterface_->getDependentVariablesSettings( )->dependentVariables_.size( ) ; i++ )
            {
                if ( itr->first == getDependentVariableId( multiArcInterface_->getDependentVariablesSettings( )
                                                           ->dependentVariables_[ i ] ) )
                {
                    sizeCurrentDependentVariable = getDependentVariableSaveSize(
                                multiArcInterface_->getDependentVariablesSettings( )->dependentVariables_[ i ] );
                }
            }

            dependentVariables.segment( currentMultiArcDependentVariableIndex, sizeCurrentDependentVariable /*multiArcDependentVariablesSize_ */) =
                    multiArcDependentVariables.segment( itr->second, sizeCurrentDependentVariable );

            currentMultiArcDependentVariableIndex += sizeCurrentDependentVariable;
        }



//        dependentVariables.block(
//                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize_, originalMultiArcStateSize_ ) =
//                multiArcDependentVariables.block(
//                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize_, originalMultiArcStateSize_ );

//        // Get single-arc matrices at current arc start
//        Eigen::MatrixXd singleArcStateTransitionAtArcStart = singleArcInterface_->getCombinedStateTransitionAndSensitivityMatrix(
//                    currentArc.second );

//        // Set coupled block
//        dependentVariables.block(
//                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ ) =
//                multiArcDependentVariables.block(
//                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ ) *
//                singleArcStateTransitionAtArcStart.block(
//                    0, 0, singleArcStateSize_, singleArcStateSize_ );

//        // Set multi-arc sensitivity block
//        dependentVariables.block(
//                    singleArcStateSize_, multiArcStateSize_, originalMultiArcStateSize_, sensitivityMatrixSize_ ) =
//                multiArcDependentVariables.block(
//                    singleArcStateSize_, multiArcStateSize_, originalMultiArcStateSize_, sensitivityMatrixSize_ );

//        std::vector< std::pair< int, int > > statePartialAdditionIndices =
//                multiArcInterface_->getStatePartialAdditionIndices( currentArc.first );

//        for( unsigned int i = 0; i < statePartialAdditionIndices.size( ); i++ )
//        {
//            combinedStateTransitionMatrix.block(
//                        statePartialAdditionIndices.at( i ).first, multiArcStateSize_,
//                        6, sensitivityMatrixSize_ ) +=
//                    combinedStateTransitionMatrix.block(
//                        statePartialAdditionIndices.at( i ).second, multiArcStateSize_,
//                        6, sensitivityMatrixSize_ );

//        }
    }

    return dependentVariables;

}

////! Function to get the concatenated state transition matrices for each arc and sensitivity matrix at a given time.
//Eigen::MatrixXd HybridArcDependentVariableInterface::getFullDependentVariable(
//        const double evaluationTime )
//{
//    Eigen::MatrixXd combinedStateTransitionMatrix = getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
//    Eigen::MatrixXd fullCombinedStateTransitionMatrix = Eigen::MatrixXd::Zero(
//                stateTransitionMatrixSize_, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_ + sensitivityMatrixSize_ );

//    fullCombinedStateTransitionMatrix.block(
//                0, 0, singleArcStateSize_, singleArcStateSize_ ) =
//            combinedStateTransitionMatrix.block(
//                0, 0, singleArcStateSize_, singleArcStateSize_ );

//    fullCombinedStateTransitionMatrix.block(
//                0, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ ) =
//            combinedStateTransitionMatrix.block(
//                0, multiArcStateSize_, singleArcStateSize_, sensitivityMatrixSize_ );

//    int currentArc = multiArcInterface_->getCurrentArc( evaluationTime ).first;

//    // Set Phi and S matrices of current arc.
//    if( currentArc >= 0 )
//    {

//        // Set multi-arc block
//        fullCombinedStateTransitionMatrix.block(
//                    singleArcStateSize_, singleArcStateSize_ + currentArc * originalMultiArcStateSize_,
//                    originalMultiArcStateSize_, originalMultiArcStateSize_ ) =
//                combinedStateTransitionMatrix.block(
//                    singleArcStateSize_, singleArcStateSize_, originalMultiArcStateSize_, originalMultiArcStateSize_ );


//        // Set coupled block
//        fullCombinedStateTransitionMatrix.block(
//                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ ) =
//                combinedStateTransitionMatrix.block(
//                    singleArcStateSize_, 0, originalMultiArcStateSize_, singleArcStateSize_ );

//        // Set multi-arc sensitivity block
//        fullCombinedStateTransitionMatrix.block(
//                    singleArcStateSize_, singleArcStateSize_ + numberOfMultiArcs_ * originalMultiArcStateSize_,
//                    originalMultiArcStateSize_, sensitivityMatrixSize_ ) =
//                combinedStateTransitionMatrix.block(
//                    singleArcStateSize_, multiArcStateSize_, originalMultiArcStateSize_, sensitivityMatrixSize_ );

//    }

//    return fullCombinedStateTransitionMatrix;
//}

}

}

