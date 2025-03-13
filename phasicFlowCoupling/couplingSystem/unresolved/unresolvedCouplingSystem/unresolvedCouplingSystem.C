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

#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::unresolvedCouplingSystem::unresolvedCouplingSystem
(
	word shapeTypeName, 
	Foam::fvMesh& mesh,
	int argc, 
	char* argv[]
)
:
	couplingSystem(shapeTypeName, mesh, argc, argv)
{
	// it searchs for cellDistribution dictionary, if not found, 
	// set the defualt cell distribtuion method 
	if( !unresolvedDict().isDict("cellDistribution") )
	{ 
		auto defDict = Foam::dictionary("cellDistribution");
		defDict.add("type", "self");
		this->subDict("unresolved").add("cellDistribution",defDict);
	}
}



pFlow::uniquePtr<pFlow::coupling::unresolvedCouplingSystem> 
	pFlow::coupling::unresolvedCouplingSystem::create
	(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]
	)
{
	
	Foam::IOdictionary couplingDict
    (
        IOobject
        (
            "couplingProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Foam::word cellDistMethod;

    if(!couplingDict.subDict("unresolved").isDict("cellDistribution"))
    {
    	cellDistMethod = "self";
    }
    else
    {
    	cellDistMethod = lookupDict<Foam::word>(couplingDict.subDict("unresolved").subDict("cellDistribution"), "type");
    }

    word couplingType = angleBracketsNames(
    	shapeTypeName +"UnresolvedCouplingSystem",
    	cellDistMethod);

    if( wordvCtorSelector_.search(couplingType))
	{
		return wordvCtorSelector_[couplingType] (shapeTypeName, mesh, argc, argv);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< couplingType << " dose not exist"
				" for drag method in "<< couplingDict.name()
				<<"\nAvaiable ones are: \n"
				,
				wordvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;

}