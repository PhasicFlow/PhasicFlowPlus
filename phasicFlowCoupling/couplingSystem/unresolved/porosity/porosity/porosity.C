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

#include "porosity.hpp"
#include "procVectorPlus.hpp"
#include "procCommunicationPlus.hpp"
#include "streams.hpp"
#include "unresolvedCouplingSystem.hpp"


pFlow::coupling::porosity::porosity
(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam
)
:
	Foam::volScalarField
	(
		Foam::IOobject
	    (
	        "alpha",
	        Foam::timeName(cMesh.mesh().time()),
	        cMesh.mesh(),
	        Foam::IOobject::MUST_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	cMesh.mesh()
	),
	alphaMin_(lookupDict<Foam::scalar>(CS.unresolvedDict().subDict("porosity"), "alphaMin")),
	uCS_(CS),
	particleDiameter_(parDiam),
	distribution_(CS.distribution())
{

}

const Foam::fvMesh &pFlow::coupling::porosity::mesh() const
{
    return uCS_.cMesh().mesh();
}


const pFlow::coupling::couplingMesh &pFlow::coupling::porosity::cMesh() const
{
    return uCS_.cMesh();
}

const pFlow::Plus::procCMField<Foam::label> &pFlow::coupling::porosity::parCellIndex() const
{
	return uCS_.parCellIndex();
}

pFlow::int32 pFlow::coupling::porosity::numInMesh() const
{
    return cMesh().numInMesh();
}

void pFlow::coupling::porosity::calculatePorosity()
{
	this->internalFieldUpdate();
	this->correctBoundaryConditions();	
}


pFlow::uniquePtr<pFlow::coupling::porosity> 
pFlow::coupling::porosity::create
(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam
)
{

	Foam::word method = 
		CS.unresolvedDict().subDict("porosity").getOrDefault<Foam::word>("method", "distribution");

	if(method == "distribution")
	{
		method = "porosity<distribution>";
	}
		

	if( dictionaryvCtorSelector_.search(method))
	{
		Foam::Info<<"    Creating porosity method "<< Green_Text(method)<<" ...\n\n";
		return dictionaryvCtorSelector_[method] (CS, cMesh, parDiam);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< method << " dose not exist"
				" for porosity method in "<< CS.unresolvedDict().subDict("porosity").name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}
