/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------  
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

/**
 * @class distributionBase
 * @brief Abstract base class for particle-to-cell distribution methods.
 *
 * Defines the interface and common functionality for all particle distribution 
 * methods used in unresolved Eulerian-Lagrangian coupling. Derived classes implement 
 * specific algorithms (PCM, Gaussian, diffusion, etc.) to map particle properties 
 * to fluid cells and vice versa.
 *
 * **Key Responsibilities:**
 * - Manage particle weights in neighboring cells
 * - Distribute particle properties to cells (particle-to-cell)
 * - Inverse-distribute cell properties to particles (cell-to-particle)
 * - Handle field smoothing operations
 * - Thread-safe value distribution with OpenMP support
 */

#ifndef __distributionBase_hpp__
#define __distributionBase_hpp__

// from std
#include <vector>

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "typeInfo.hpp"
#include "virtualConstructor.hpp"

// from phasicFlowPlus
#include "procCMField.hpp"

namespace pFlow::coupling
{

class couplingMesh;


class distributionBase
{
public:
    
    /// Type alias for cell index and weight pair
    using cellWeight = std::pair<Foam::label, Foam::scalar>;

private:
    
    /// Flag indicating whether cell distribution is used
    const bool useCelldistribution_;

    /// Particle weights in neighboring cells
    Plus::procCMField<std::vector<cellWeight>>          weights_;

    /// Reference to the coupling mesh
    const couplingMesh&	                                cMesh_;

protected:
    
    /// Get mutable reference to particle weights
    inline 
    Plus::procCMField<std::vector<cellWeight>>& weights()
    {
        return weights_;
    }

    /// Get const reference to all particle weights
    inline 
    const Plus::procCMField<std::vector<cellWeight>>& weights()const
    {
        return weights_;
    }

    /// Get mutable reference to weights for a specific particle
    inline 
    std::vector<cellWeight>& weights(Foam::label parIndx)
    {
        return weights_[parIndx];
    }

    /// Get const reference to weights for a specific particle
    inline 
    const std::vector<cellWeight>& weights(Foam::label parIndx)const
    {
        return weights_[parIndx];
    }

    /// Construct neighbor lists with search length and max layers (pure virtual)
    virtual 
    void constructLists(
        const Foam::scalar searchLen, 
        const Foam::label maxLayers) = 0;
    
    /// Construct neighbor lists with max layers only (pure virtual)
    virtual 
    void constructLists(const Foam::label maxLayers)=0;

public:
    
    /// Type info
    TypeInfo("distributionBase");

    /// Constructor with dictionary for initialization
    distributionBase(
        bool                         useCellDistribution,
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Constructor without dictionary
    distributionBase(
        bool                         useCellDistribution, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~distributionBase() = default;

    create_vCtor
    (
        distributionBase,
        dictionary,
        (
            const Foam::dictionary& 	 parrentDict, 
            const couplingMesh& 	     cMesh,
            const Plus::centerMassField& centerMass
        ),
        (parrentDict, cMesh, centerMass)
    );

    /// Check if cell distribution is required
    inline bool requireCellDistribution()const
    {
        return useCelldistribution_;
    }

    /// Get reference to the coupling mesh
    inline
    const couplingMesh& cMesh()const
    {
        return cMesh_;
    }

    /// Get reference to the underlying OpenFOAM mesh
    const Foam::fvMesh& mesh()const;
    
    /// Distribute scalar value to cells (thread-safe with OpenMP)
    inline
    void distributeValue_OMP(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volScalarField::Internal& internalField,
        const Foam::scalar& val)const
    {
        if(useCelldistribution_)
        {
            const auto& weightsPar = weights_[parIndx];
            for(const auto& [cellIndx, weight]: weightsPar)
            {
                #pragma omp atomic
                internalField[cellIndx] += val*weight;
            }
        }
        else
        {
            #pragma omp atomic
            internalField[parCellIndx] += val;
        } 
    }

    /// Distribute vector value to cells (thread-safe with OpenMP)
    inline
    void distributeValue_OMP(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volVectorField::Internal& internalField,
        const Foam::vector& val )const
    {
        if(useCelldistribution_)
        {
            const auto& weightsPar = weights_[parIndx];        
            for(const auto& [cellIndx, weight]: weightsPar)
            {
                const auto v = weight* val;
                auto& tv = internalField[cellIndx]; 
                #pragma omp atomic
                tv.x() += v.x();
                
                #pragma omp atomic
                tv.y() += v.y();

                #pragma omp atomic
                tv.z() += v.z();
            }
        }
        else
        {
            auto& tv = internalField[parCellIndx]; 
            #pragma omp atomic
            tv.x() += val.x();
            
            #pragma omp atomic
            tv.y() += val.y();

            #pragma omp atomic
            tv.z() += val.z();
        }
    }

    /// Distribute scalar value to cells (non-threaded)
    inline 
    void distributeValue(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volScalarField::Internal& internalField,
        const Foam::scalar& val)const
    {
        if(useCelldistribution_)
        {
            const auto& weightsPar = weights_[parIndx];
            for(const auto& [cellIndx, weight]: weightsPar)
            {
                internalField[cellIndx] += val*weight;
            }
        }
        else
        {
            internalField[parCellIndx] += val;
        }
    }

    /// Distribute vector value to cells (non-threaded)
    inline
    void distributeValue(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volVectorField::Internal& internalField,
        const Foam::vector& val)const
    {
        if(useCelldistribution_)
        {
            const auto& weightsPar = weights_[parIndx];
                    
            for(const auto& [cellIndx, weight]: weightsPar)
            {
                internalField[cellIndx] += weight* val;
            }
        }
        else
        {
            internalField[parCellIndx] += val;
        }
    }

    /// Inverse-distribute vector field value to particle (cell-to-particle)
    inline 
    void inverseDistributeValue(
        Foam::label parIndx,
        Foam::label parCellIndx,
        const Foam::volVectorField::Internal& internalField,
        Foam::vector& val)const
    {
        if(useCelldistribution_)
        {
            Foam::vector avVal =Foam::vector(0,0,0);
            const auto& weightsPar = weights_[parIndx];
            for(const auto& [cellIndx, w]: weightsPar)
            {
                avVal += w*internalField[cellIndx];
            }
            val = avVal;
        }
        else
        {
            val = internalField[parCellIndx];
        }
    }

    /// Inverse-distribute scalar field value to particle (cell-to-particle)
    inline 
    void inverseDistributeValue(
        Foam::label parIndx,
        Foam::label parCellIndx,
        const Foam::volScalarField::Internal& internalField,
        Foam::scalar& val)const
    {
        if(useCelldistribution_)
        {
            Foam::scalar avVal = 0.0;

            const auto& weightsPar = weights_[parIndx];
            for(const auto& [cellIndx, w]: weightsPar)
            {
                avVal += w*internalField[cellIndx];
            }
            val = avVal;
        }
        else
        {
            val = internalField[parCellIndx];
        }
    }
    
    /// Update distribution weights for all particles (pure virtual)
    virtual 
    void updateWeights(const Plus::procCMField<real> & parDiameter) = 0;

    /// Smooth a vector field using the distribution method (pure virtual)
    virtual 
    void smoothenField(Foam::volVectorField& field)const = 0;
    
    /// Smooth a scalar field using the distribution method (pure virtual)
    virtual
    void smoothenField(Foam::volScalarField& field)const = 0;

    /// Get the name of the distribution method (pure virtual)
    virtual 
    Foam::word distributionMethodName()const = 0;

    /// Factory method to create a distribution instance from dictionary
    static
    uniquePtr<distributionBase> create
    (
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass
    );
	    
}; // end distributionBase


} // end pFlow::coupling

#endif //__distributionBase_hpp__
