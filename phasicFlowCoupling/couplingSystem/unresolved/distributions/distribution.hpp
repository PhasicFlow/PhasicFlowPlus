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
#ifndef __distribution_hpp__
#define __distribution_hpp__

// from std
#include <set>
#include <vector>

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "typeInfo.hpp"

// from phasicFlowPlus
#include "procCMField.hpp"

namespace pFlow::coupling
{

class couplingMesh;


class distribution
{
public:
    
    using cellWeight = std::pair<Foam::label, Foam::scalar>;

protected:

    Plus::procCMField<std::vector<cellWeight>>          weights_;

    // list of neighbor cells for each cell 
    std::vector<std::vector<Foam::label>>               neighborList_;
    
    const Foam::fvMesh&	                                mesh_;

    // members

    void constructLists(const Foam::scalar searchLen, const Foam::label maxLayers=3);
    
    void constructLists(const Foam::label maxLayers = 3);

    void parseNeighbors(
		const Foam::label 		targetCelli,
		const Foam::vector& 	targetCellCentre, 
		const Foam::scalar 		searchLen,
		const Foam::scalar 		celli,
		std::set<Foam::label>& 	finalList,
		const Foam::label 		layerNumber,
        const Foam::label 		maxLayers=3);

    void parseNeighbors(
        const Foam::label       targetCelli,
        const Foam::label       layerNumber,
        const Foam::label       maxLayers, 
        const Foam::label       celli, 
        std::set<Foam::label>&  finalList       
    );

public:
    
/// Type info
    TypeInfoNV("distribution");

    /// Construct from dictionary 
    distribution(
        const Foam::dictionary& 	 dict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    ~distribution() = default;

    inline
    void distributeValue_OMP(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volScalarField::Internal& internalField,
        const Foam::scalar& val)const
    {
        const auto& weightsPar = weights_[parIndx];
        for(const auto& [cellIndx, weight]: weightsPar)
        {
            #pragma omp atomic
            internalField[cellIndx] += val*weight;
        }
    }

    inline
    void distributeValue_OMP(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volVectorField::Internal& internalField,
        const Foam::vector& val )const
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

    inline 
    void distributeValue(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volScalarField::Internal& internalField,
        const Foam::scalar& val)const
    {
        const auto& weightsPar = weights_[parIndx];
        for(const auto& [cellIndx, weight]: weightsPar)
        {
            internalField[cellIndx] += val*weight;
        }
    }

    inline
    void distributeValue(
        Foam::label parIndx, 
        Foam::label parCellIndx,
        Foam::volVectorField::Internal& internalField,
        const Foam::vector& val)const
    {
        const auto& weightsPar = weights_[parIndx];
        for(const auto& [cellIndx, weight]: weightsPar)
        {
            internalField[cellIndx] += val*weight;
        }
    }

    inline 
    void inverseDistributeValue(
        Foam::label parIndx,
        Foam::label parCellIndx,
        const Foam::volVectorField::Internal& internalField,
        Foam::vector& val 
    )const
    {
        Foam::vector avVal =Foam::vector(0,0,0);

        const auto& weightsPar = weights_[parIndx];
        for(const auto& [cellIndx, w]: weightsPar)
        {
            avVal += w*internalField[cellIndx];
        }
        
        val = avVal;
    }

    inline 
    void inverseDistributeValue(
        Foam::label parIndx,
        Foam::label parCellIndx,
        const Foam::volScalarField::Internal& internalField,
        Foam::scalar& val 
    )const
    {
        Foam::scalar avVal = 0.0;

        const auto& weightsPar = weights_[parIndx];
        for(const auto& [cellIndx, w]: weightsPar)
        {
            avVal += w*internalField[cellIndx];
        }
        
        val = avVal;
    }

    void smoothenField(Foam::volVectorField& field)const
    {}

    void smoothenField(Foam::volScalarField& field)const
	{}

    const Plus::procCMField<std::vector<cellWeight>> weights()const
    {
        return weights_;
    }
    
}; // end distribution


} // end pFlow::coupling

#endif //__distribution_hpp__
