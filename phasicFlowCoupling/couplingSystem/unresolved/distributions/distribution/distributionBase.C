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

#include "distributionBase.hpp"
#include "processorPlus.hpp"
#include "couplingMesh.hpp"

pFlow::coupling::distributionBase::distributionBase(
    bool                         useCellDistribution,
    const Foam::dictionary& 	 parrentDict, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    useCelldistribution_(useCellDistribution),
    weights_("weights", centerMass),
    cMesh_(cMesh)
{}

pFlow::coupling::distributionBase::distributionBase(
    bool                         useCellDistribution, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    useCelldistribution_(useCellDistribution),
    weights_("weights", centerMass),
    cMesh_(cMesh)
{}

const Foam::fvMesh &pFlow::coupling::distributionBase::mesh() const
{
    return cMesh_.mesh();
}

pFlow::uniquePtr<pFlow::coupling::distributionBase> 
    pFlow::coupling::distributionBase::create
(
    const Foam::dictionary &parrentDict, 
    const couplingMesh &cMesh, 
    const Plus::centerMassField &centerMass
)
{
    auto distType = parrentDict.get<Foam::word>("distributionMethod");

    if( dictionaryvCtorSelector_.search(distType))
    {
        Foam::Info<<"    Crearing cell distribution "<<Green_Text(distType)<<" for coupling system...\n\n";
        return dictionaryvCtorSelector_[distType]
        (
            parrentDict,
            cMesh,
            centerMass
        );
    }
    else
    {
        if(Plus::processor::isMaster())
        {
            printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< distType << " dose not exist"
				" for distribution method in "<< parrentDict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
        }
        Plus::processor::abort(0);
    }

    return nullptr;
}
