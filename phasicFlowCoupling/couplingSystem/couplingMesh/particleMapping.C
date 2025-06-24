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

#include "particleMapping.hpp"
#include "couplingMesh.hpp"
#include "procDEMSystemPlus.hpp"


pFlow::coupling::particleMapping::particleMapping
(
    const Foam::dictionary &dict
)
:
    Plus::procCommunication(),
    centerMass_(),
    domainExpansionRatio_
    (
        Foam::max(lookupDict<Foam::scalar>(dict, "domainExpansionRatio"), 0.5)
    ),
    domainUpdateInterval_
    (
        lookupDict<Foam::scalar>(dict, "domainUpdateInterval")
    )
{

}


bool pFlow::coupling::particleMapping::checkForDomainUpdate
(
    Foam::scalar t, 
    Foam::scalar fluidDt
)const
{
    if( !firstConstructed_ )
    {
        return true;
    }

    if( std::abs(t-lastTimeUpdated_) < static_cast<Foam::scalar>(0.98*fluidDt) )
    {
        return true;
    }
    
    if( std::abs(t-(lastTimeUpdated_+domainUpdateInterval_)) < static_cast<Foam::scalar>(0.98*fluidDt))
    {
        return true;
    }

    return false;
}

bool pFlow::coupling::particleMapping::update
(
    Foam::scalar t, 
    Foam::scalar fluidDt, 
    Plus::procDEMSystem &pDEMSystem, 
    const couplingMesh &cMesh
)
{
    
    if( checkForDomainUpdate(t, fluidDt) )
    {
        firstConstructed_ = true;

        lastTimeUpdated_ = t;

        REPORT(0)<<Blue_Text("Particle mapping in processors at time :")<< 
            Yellow_Text(t) <<" s"<<END_REPORT;

        auto mBox = cMesh.meshBox();

        if(! collectAllToAll(mBox, meshBoxes_))
        {
            fatalErrorInFunction<<"could not distribute meshBox"<<endl;
            Plus::processor::abort(0);
        }

        REPORT(1)<< "Mesh boxes on all processors updated"<<pFlow::endl;

        if(!pDEMSystem.updateParticleDistribution(
            domainExpansionRatio_, 
            meshBoxes_))
        {
            fatalErrorInFunction;
            Plus::processor::abort(0);
            return false;
        }

        REPORT(1)<< "Re-mapped particles on all boxes"<<pFlow::endl;

        auto numParsInDomains = pDEMSystem.numParInDomainMaster();

        if( auto [thisNoPars, success] =  
            distributeMasterToAll(numParsInDomains); !success)
        {
            fatalErrorInFunction<<
            "failed to distribute particle numbers among processors"<<endl;
            Plus::processor::abort(0);
            return false;
        }
        else
        {
            if(!centerMass_.checkForNewSize(thisNoPars))
            {
                fatalErrorInFunction<<
                "cannot change the size of containers to new size "<< 
                thisNoPars<<endl;
                Plus::processor::abort(0);
                return false;
            }
        }

        // first cunstructs index distribution
        auto parIndexInDomains = pDEMSystem.parIndexInDomainsMaster();

        if(!realScatteredComm_.changeDataMaps(parIndexInDomains))
        {
            fatalErrorInFunction<<
            "error in creating index block for real type"<<endl;
            Plus::processor::abort(0);
            return false;
        }

        if(!realx3ScatteredComm_.changeDataMaps(parIndexInDomains))
        {
            fatalErrorInFunction<<
            "error in creating index block for realx3 type"<<endl;
            Plus::processor::abort(0);
            return false;
        }

        if(!uint32ScatteredComm_.changeDataMaps(parIndexInDomains))
        {
            fatalErrorInFunction<<
            "error in creating index block for realx3 type"<<endl;
            Plus::processor::abort(0);
            return false;
        }

        REPORT(1)<< "Data mapping updated"<<pFlow::endl;
    }

    auto allPos = pDEMSystem.particlesCenterMassAllMaster();
	auto thisPos = makeSpan(centerMass_);
	if(!realx3ScatteredComm_.distribute(allPos, thisPos))
	{
		fatalErrorInFunction<<
		"cannot distribute particle positions to processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}

    return true;
}
