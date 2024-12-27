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
#include "Timer.hpp"


void pFlow::coupling::porosity::mapCenters()
{
	numInMesh_ = 0;
	size_t numPar = centerMass_.size();

	
	#pragma omp parallel for reduction(+:numInMesh_)
	for(size_t i = 0; i<numPar; i++)
	{
		auto cellId = cMesh_.findCellTree(centerMass_[i], parCellIndex_[i]);
		parCellIndex_[i] = cellId;
		if( cellId >= 0 ) numInMesh_++;	
	}
}


pFlow::coupling::porosity::porosity(
	Foam::dictionary		dict, 
	couplingMesh& 			cMesh, 
	Plus::centerMassField& 	centerMass, 
	Plus::realProcCMField& 	parDiam)
:
	Foam::volScalarField
	(
		Foam::IOobject
	    (
	        "alpha",
	        cMesh.mesh().time().timeName(),
	        cMesh.mesh(),
	        Foam::IOobject::MUST_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	cMesh.mesh()
	),
	alphaMin_(dict.lookup<Foam::scalar>("alphaMin")),
	cMesh_(cMesh),
	centerMass_(centerMass),
	particleDiameter_(parDiam),
	parCellIndex_
	(
		"parCellIndex",
		static_cast<Foam::label>(-1),
		centerMass,
		true
	)
{

}

void pFlow::coupling::porosity::calculatePorosity()
{
	Timer t;
	t.start();
	mapCenters();
	this->internalFieldUpdate();
	t.end();
	//output<<"mapping execution time " << t.lastTime()<<endl;
	
	this->correctBoundaryConditions();
	
}

void pFlow::coupling::porosity::reportNumInMesh()
{
	Plus::procCommunication comm;
	if( auto [numInMeshAll, success] = comm.collectAllToMaster(numInMesh()); success)
	{
		if(Plus::processor::isMaster())
		{
			int32 s=0;
			for(auto v:numInMeshAll) s += v;

			output<<Blue_Text("Particles located in processor meshes:") << Yellow_Text(numInMeshAll)<<
			" => "<< Yellow_Text(s)<< endl;
		}
	}
}

pFlow::uniquePtr<pFlow::coupling::porosity> 
pFlow::coupling::porosity::create
(
	Foam::dictionary		dict, 
	couplingMesh& 			cMesh, 
	Plus::centerMassField& 	centerMass, 
	Plus::realProcCMField& 	parDiam
)
{
	auto method = dict.lookup<Foam::word>("method");
	if( dictionaryvCtorSelector_.search(method))
	{
		return dictionaryvCtorSelector_[method] (dict, cMesh, centerMass, parDiam);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< method << " dose not exist"
				" for porosity method in "<< dict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}
