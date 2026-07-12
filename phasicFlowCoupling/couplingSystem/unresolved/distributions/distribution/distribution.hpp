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
#include "distributionBase.hpp"

namespace pFlow::coupling
{

class distribution
:
    public distributionBase
{    
protected:

    // list of neighbor cells for each cell 
    std::vector<std::vector<Foam::label>>               neighborList_;
    
    // members functions
    void parseNeighbors(
		const Foam::label 		targetCelli,
		const Foam::vector& 	targetCellCentre, 
		const Foam::scalar 		searchLen,
		const Foam::scalar 		celli,
		std::set<Foam::label>& 	finalList,
		const Foam::label 		layerNumber,
        const Foam::label 		maxLayers);

    void parseNeighbors(
        const Foam::label       targetCelli,
        const Foam::label       layerNumber,
        const Foam::label       maxLayers, 
        const Foam::label       celli, 
        std::set<Foam::label>&  finalList);
    
    void constructLists(const Foam::scalar searchLen, const Foam::label maxLayers)override;
    
    void constructLists(const Foam::label maxLayers)override;


public:
    
/// Type info
    TypeInfo("distribution");

    /// Construct from dictionary 
    distribution(
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~distribution() = default;

    void smoothenField(Foam::volScalarField& field)const override;

    void smoothenField(Foam::volVectorField& field)const override;
    
}; // end distribution


} // end pFlow::coupling

#endif //__distribution_hpp__
