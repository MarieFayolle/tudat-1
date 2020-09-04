/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DEPENDENTVARIABLESINTERFACE_H
#define TUDAT_DEPENDENTVARIABLESINTERFACE_H

#include <iostream>
#include <vector>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutput.h"

namespace tudat
{

namespace propagators
{

//! Base class for interface object of interpolation of numerically propagated dependent variables.
/*!
 *  Base class for interface object of interpolation of numerically propagated dependent variables.
 *  Derived classes implement the case of single-arc/multi-arc/hybrid combined dynamics.
 */
class DependentVariablesInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * state transition matrix.
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     */
    DependentVariablesInterface(
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesSettings/*,
            const int numberOfInitialDynamicalParameters,
            const int numberOfParameters*/ ):
    dependentVariablesSettings_( dependentVariablesSettings )
    {
        singleDependentVariableSettingsVector_ = dependentVariablesSettings_->dependentVariables_;
        dependentVariablesSize_ = 0;

        dependentVariablesIdsAndIndices_.clear( );

        for ( unsigned int i= 0 ; i < singleDependentVariableSettingsVector_.size( ) ; i++ )
        {
            dependentVariablesTypes_.push_back( singleDependentVariableSettingsVector_[ i ]
                                                ->dependentVariableType_ );

            dependentVariablesIdsAndIndices_[ getDependentVariableId( singleDependentVariableSettingsVector_[ i ] ) ]
                    = dependentVariablesSize_;

            dependentVariablesSize_ += getDependentVariableSaveSize( singleDependentVariableSettingsVector_[ i ] );

        }
    }

    //! Destructor.
    virtual ~DependentVariablesInterface( ){ }

    //! Function to get the dependent variables at a given time.
    /*!
     *  Function to get the dependent variabless at a given time.
     *  \param evaluationTime Time at which to evaluate dependent variables interpolator
     *  \return Concatenated dependent variables.
     */
    virtual Eigen::VectorXd getDependentVariables( const double evaluationTime ) = 0;

    //! Function to get the value of a single dependent variable at a given time.
    Eigen::VectorXd getSingleDependentVariable(
            const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( getDependentVariableSaveSize( dependentVariableSettings ) );

        // Retrieve ID and index of dependent variable of interest.
        std::string dependentVariableId = getDependentVariableId( dependentVariableSettings );
//        std::map< std::string, int >::iterator iterator = dependentVariablesIdsAndIndices_.find( dependentVariableId );
        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

        // Retrieve full vector of dependent variables at a given time.
        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

        dependentVariable =
                fullDependentVariablesVector.segment( dependentVariableIndex, getDependentVariableSaveSize( dependentVariableSettings ) );

        return dependentVariable;
    }

    //! Function to get the value of a single dependent variable at a given time, from the dependent variable ID.
    Eigen::VectorXd getSingleDependentVariable(
            const std::string dependentVariableId,
            const int dependentVariableSize,
            const double evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( dependentVariableSize );

        // Retrieve index of dependent variable of interest.
        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

        // Retrieve full vector of dependent variables at a given time.
        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

        dependentVariable =
                fullDependentVariablesVector.segment( dependentVariableIndex, dependentVariableSize );

        return dependentVariable;
    }



    //! Function to get the size of the dependent variables
    /*!
     * Function to get the size of the dependent variables
     * \return Size of dependent variables
     */
    int getDependentVariablesize( )
    {
        return dependentVariablesSize_;
    }

    //! Function to retrieve the dependent variable settings object.
    /*!
     * Function to retrieve the dependent variable settings object.
     * \return dependent variable settings
     */
    std::shared_ptr< propagators::DependentVariableSaveSettings > getDependentVariablesSettings( )
    {
        return dependentVariablesSettings_;
    }

    //! Function to retrieve vector of single dependent variable settings.
    /*!
     * Function to retrieve the dependent variable settings object.
     * \return dependent variable settings
     */
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > >
    getVectorSingleDependentVariableSettings( )
    {
        return singleDependentVariableSettingsVector_;
    }


   //! Function to retrieve the map containing the dependent variables Ids and indices.
   std::map< std::string, int > getDependentVariablesIdsAndIndices( )
   {
       return dependentVariablesIdsAndIndices_;
   }


protected:

    //! Dependent variable settings object
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesSettings_;

    //! Vector of single dependent variable settings objects
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > singleDependentVariableSettingsVector_;

    //! Type of the dependent variables of interest
    std::vector< propagators::PropagationDependentVariables > dependentVariablesTypes_;

    //! Size of dependent variable vector
    int dependentVariablesSize_;

    //! Map containing the dependent variables Ids and indices
    std::map< std::string, int > dependentVariablesIdsAndIndices_;

};

//! Interface object of interpolation of numerically propagated dependent variable for single-arc
//! estimation.
class SingleArcDependentVariablesInterface : public DependentVariablesInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateTransitionMatrixInterpolator Interpolator returning the state transition matrix as a function of time.
     * \param sensitivityMatrixInterpolator Interpolator returning the sensitivity matrix as a function of time.
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * state transition matrix.
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    SingleArcDependentVariablesInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >
            dependentVariablesInterpolator,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesSettings ):
        DependentVariablesInterface( dependentVariablesSettings ),
        dependentVariablesInterpolator_( dependentVariablesInterpolator )
    {
        dependentVariables_ = Eigen::VectorXd::Zero( dependentVariablesSize_ );
    }

    //! Destructor.
    ~SingleArcDependentVariablesInterface( ){ }

    //! Function to reset the dependent variable interpolators
    /*!
     * Function to reset the dependent variable interpolators
     * \param stateTransitionMatrixInterpolator New interpolator returning the state transition matrix as a function of time.
     * \param sensitivityMatrixInterpolator New interpolator returning the sensitivity matrix as a function of time.
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    void updateDependentVariablesInterpolators(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >
            dependentVariablesInterpolator );

    //! Function to get the interpolator returning the dependent variable as a function of time.
    /*!
     * Function to get the interpolator returning the dependent variable as a function of time.
     * \return Interpolator returning the dependent variable as a function of time.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >
    getDependentVariablesInterpolator( )
    {
        return dependentVariablesInterpolator_;
    }

    //! Function to get the concatenated state transition and sensitivity matrix at a given time.
    /*!
     *  Function to get the concatenated state transition and sensitivity matrix at a given time.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Concatenated state transition and sensitivity matrices.
     */
    Eigen::VectorXd getDependentVariables( const double evaluationTime );


private:

    //! Predefined vector to use as return value when calling getDependentVariables.
    Eigen::VectorXd dependentVariables_;

    //! Interpolator returning the dependent variables as a function of time.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >
    dependentVariablesInterpolator_;


};

//! Interface object of interpolation of numerically propagated dependent variable for multi-arc
//! estimation.
class MultiArcDependentVariablesInterface: public DependentVariablesInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateTransitionMatrixInterpolators Interpolators returning the state transition matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param sensitivityMatrixInterpolators Interpolators returning the sensitivity matrix as a function of time, vector
     * entries represent matrix history for each arc.
     * \param arcStartTimes Times at which the multiple arcs start
     * \param arcEndTimes Times at which the multiple arcs end
     * \param numberOfInitialDynamicalParameters Size of the estimated initial state vector (and size of square
     * sing-arc state transition matrix times number of arcs.)
     * \param numberOfParameters Total number of estimated parameters (initial states and other parameters).
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    MultiArcDependentVariablesInterface(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
            dependentVariablesInterpolators,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesSettings,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes );

    //! Destructor
    ~MultiArcDependentVariablesInterface( ){ }

    //! Function to reset the dependent variable interpolators
    /*!
     * Function to reset the dependent variable interpolators
     * \param dependentVariableInterpolators New vector of interpolators returning the dependent variable as a
     * function of time.
     * \param arcStartTimes Times at which the multiple arcs start.
     * \param arcEndTimes Times at which the multiple arcs end
     * \param statePartialAdditionIndices Vector of pair providing indices of column blocks of variational equations to add to other column blocks
     */
    void updateDependentVariablesInterpolators(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
            dependentVariablesInterpolators,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes  );

    //! Function to get the vector of interpolators returning the dependent variable as a function of time.
    /*!
     * Function to get the vector of interpolators returning the dependent variable as a function of time.
     * \return Vector of interpolators returning the dependent variable as a function of time.
     */
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
    getDependentVariablesInterpolators( )
    {
        return dependentVariablesInterpolators_;
    }

    //! Function to get the concatenated single-arc dependent variable at a given time.
    /*!
     *  Function to get the concatenated single-arc dependent variable at a given time, evaluates matrices
     *  at the arc in which evaluationTime is located.
     *  \param evaluationTime Time at which to evaluate dependent variable interpolators
     *  \return Concatenated dependent variable.
     */
    Eigen::VectorXd getDependentVariables( const double evaluationTime );

    //! Function to retrieve the current arc for a given time
    /*!
     * Function to retrieve the current arc for a given time
     * \param evaluationTime Time at which current arc is to be determined
     * \return Pair with current arc index and associated arc initial time.
     */
    std::pair< int, double >  getCurrentArc( const double evaluationTime );

    //! Function to retrieve the number of arcs in dynamics
    /*!
     * Function to retrieve the number of arcs in dynamics
     * \return Number of arcs in dynamics
     */
    int getNumberOfArcs( )
    {
        return arcStartTimes_.size( );
    }


private:

    //! List of interpolators returning the dependent variable as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
    dependentVariablesInterpolators_;

    //! Times at which the multiple arcs start
    std::vector< double > arcStartTimes_;

    std::vector< double > arcEndTimes_;

    //! Number of arcs.
    int numberOfStateArcs_;

    //! Look-up algorithm to determine the arc of a given time.
    std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;

};

//! Interface object of interpolation of numerically propagated dependent variable for a hybrid of
//! single-and multi-arc estimation (single order is put first in concatenation)
/*!
 *  Interface object of interpolation of numerically propagated dependent variable for a hybrid of
 *  single-and multi-arc estimation (single order is put first in concatenation). The single- and multi-arc given as input must
 *  be consistent: the single arc bodies/states must also be included in the multi-arc model, in order to properly generate
 *  the coupling terms between single and multi-arc states (see HybridArcVariationalEquationsSolver).
 */
class HybridArcDependentVariablesInterface: public DependentVariablesInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param singleArcInterface Object to retrieve the dependent variable for single arc component
     * \param multiArcInterface Object to retrieve the dependent variable for multi arc component
     */
    HybridArcDependentVariablesInterface(
            const std::shared_ptr< SingleArcDependentVariablesInterface > singleArcInterface,
            const std::shared_ptr< MultiArcDependentVariablesInterface > multiArcInterface ):
        DependentVariablesInterface( multiArcInterface->getDependentVariablesSettings( ) ),
        singleArcInterface_( singleArcInterface ), multiArcInterface_( multiArcInterface )
    {
        singleArcDependentVariablesSize_ = 0;
        multiArcDependentVariablesSize_ = 0;
        numberOfMultiArcs_ = 0;

        if ( singleArcInterface_ != nullptr )
        {
            singleArcDependentVariablesIdsAndIndices_ = singleArcInterface_->getDependentVariablesIdsAndIndices( );
            singleArcDependentVariablesSize_ = singleArcInterface_->getDependentVariablesize( );
        }
        std::map< std::string, int > multiArcDependentVariablesIdsAndIndices;
        if ( multiArcInterface_ != nullptr )
        {
            multiArcDependentVariablesIdsAndIndices = multiArcInterface_->getDependentVariablesIdsAndIndices( );
            multiArcDependentVariablesSize_ = multiArcInterface_->getDependentVariablesize( ) - singleArcDependentVariablesSize_;
            numberOfMultiArcs_ = multiArcInterface->getNumberOfArcs( );
        }

        // Check input consistency: verify that the same dependent variable is not saved twice (in both the single and multi-arc
        // dependent variables interface)
        if ( ( singleArcInterface != nullptr ) && ( multiArcInterface != nullptr ) )
        {
            for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices.begin( ) ;
                  itr != multiArcDependentVariablesIdsAndIndices.end( ) ; itr++ )
            {
                if ( singleArcDependentVariablesIdsAndIndices_.count( itr->first ) != 0 )
                {

                    // Remove dependent variable from multi-arc dependent variables interface.
                    IdsAndIndicesMultiArcDependentVariablesToBeRemoved_[ itr->first ] = itr->second;

                    std::cerr << "Warning when making hybrid dependent variables interface, dependent variable "
                              << itr->first << " required by the multi-arc interface is already accounted for "
                                               "in the single arc interface." << std::endl;

//                    throw std::runtime_error( "Warning when making hybrid dependent variables interface, dependent variable " +
//                                              itr->first + "required by the multi-arc interface is already accounted for "
//                                                               "in the single arc interface." );
                }
            }
        }

//        std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > >
//                newMultiArcDependentVariablesInterpolators;
//        std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > newMultiArcDependentVariablesSettings;
        int sizeRemovedDependentVariables = 0;
        int counterDependentVariable = 0;
        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices.begin( ) ;
              itr != multiArcDependentVariablesIdsAndIndices.end( ) ; itr++ )
        {
            if ( IdsAndIndicesMultiArcDependentVariablesToBeRemoved_.count( itr->first ) == 0 )
            {
                multiArcDependentVariablesIdsAndIndices_[ itr->first ] = itr->second - sizeRemovedDependentVariables;
            }
            else
            {
                sizeRemovedDependentVariables += getDependentVariableSaveSize(
                            multiArcInterface_->getDependentVariablesSettings( )->dependentVariables_[ counterDependentVariable ] );
            }
            counterDependentVariable += 1;
//            counter += 1;
        }



//        singleArcDependentVariablesSize_ = singleArcInterface_->getDependentVariablesize( );
//        multiArcDependentVariablesSize_ = multiArcInterface_->getDependentVariablesize( );
        hybridArcDependentVariablesSize_ = singleArcDependentVariablesSize_ + multiArcDependentVariablesSize_;
//        originalMultiArcDependentVariablesSize_ = multiArcDependentVariablesSize_ - singleArcDependentVariablesSize_;

//        numberOfMultiArcs_ = multiArcInterface->getNumberOfArcs( );


        // Reset dependent variables IDs and indices map.
        /// TO BE MODIFIED AFTER REMOVING PART OF THE MULTI-ARC DEPENDENT VARIABLES IF NEEDED
        dependentVariablesIdsAndIndices_.clear( );
        dependentVariablesIdsAndIndices_ = singleArcDependentVariablesIdsAndIndices_;
        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices_.begin( ) ;
              itr != multiArcDependentVariablesIdsAndIndices_.end( ) ; itr++ )
        {
            dependentVariablesIdsAndIndices_[ itr->first ] = itr->second + singleArcDependentVariablesSize_;
        }

////        if ( multiArcInterface->getVectorSingleDependentVariableSettings( ).size( ) !=
////             singleArcInterface->getVectorSingleDependentVariableSettings( ).size( ) )
////        {
////            throw std::runtime_error( "Error when making hybrid dependent variable interface, number of dependent "
////                                      "variables is inconsistent between single arc and multi arc interfaces" );
////        }

////        for ( unsigned int i = 0 ; i < multiArcInterface->getVectorSingleDependentVariableSettings( ).size( ) ; i++ )
////        {
////            if ( getDependentVariableId( multiArcInterface->getVectorSingleDependentVariableSettings( )[ i ] )
////                 != getDependentVariableId( singleArcInterface->getVectorSingleDependentVariableSettings( )[ i ] ) )
////            {
////                throw std::runtime_error( "Error when making hybrid dependent variable interface, "
////                                          "dependent variable IDs are inconsistent between single arc "
////                                          "and multi arc interfaces." );
////            }
////        }

    }

    //! Destructor
    ~HybridArcDependentVariablesInterface( ){ }

    //! Function to get the dependent variable at a given time.
    /*!
     *  Function to get the dependent variable matrix at a given time.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Dependent variable.
     */
    Eigen::VectorXd getDependentVariables( const double evaluationTime );

    //! Function to get the value of a single dependent variable at a given time.
    Eigen::VectorXd getSingleDependentVariable(
            const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( getDependentVariableSaveSize( dependentVariableSettings ) );

        // Retrieve ID and index of dependent variable of interest.
        std::string dependentVariableId = getDependentVariableId( dependentVariableSettings );
//        std::map< std::string, int >::iterator iterator = dependentVariablesIdsAndIndices_.find( dependentVariableId );
        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

        // Retrieve full vector of dependent variables at a given time.
        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

        dependentVariable =
                fullDependentVariablesVector.segment( dependentVariableIndex, getDependentVariableSaveSize( dependentVariableSettings ) );

        return dependentVariable;
    }

private:

    //! Object to retrieve dependent variable for single arc component
    std::shared_ptr< SingleArcDependentVariablesInterface > singleArcInterface_;

    //! Object to retrieve dependent variable for multi arc component
    std::shared_ptr< MultiArcDependentVariablesInterface > multiArcInterface_;

    //! Dependent variables IDs and indices for single-arc interface.
    std::map< std::string, int > singleArcDependentVariablesIdsAndIndices_;

    //! Dependent variables IDs and indices for multi-arc interface.
    std::map< std::string, int > multiArcDependentVariablesIdsAndIndices_;

    //! IDs and indices of dependent variables to be removed from multi-arc interface (in case already included in single arc
    //! interface)
    std::map< std::string, int > IdsAndIndicesMultiArcDependentVariablesToBeRemoved_;

    //! Size of single-arc dependent variables
    int singleArcDependentVariablesSize_;

    //! Size of multi-arc dependent variables.
    int multiArcDependentVariablesSize_;

    //! Size of hybrid arc dependent variables.
    int hybridArcDependentVariablesSize_;

    //!

//    //! Full dependent variables size of single arc in original multi-arc model (full multi-arc size minus single-arc size).
//    int originalMultiArcDependentVariablesSize_;

    //! Number of arcs in multi-arc model.
    int numberOfMultiArcs_;
};



} // namespace propagators

} // namespace tudat

#endif // TUDAT_DEPENDENTVARIABLESINTERFACE_H
