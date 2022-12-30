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

#include "couplingSystem.hpp"

pFlow::coupling::couplingSystem::couplingSystem(
		word demSystemName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[])
:
	couplingMesh_(mesh),
	processorComm_(),
	procDEMSystem_(demSystemName, argc, argv),
	centerMass_(),
	particleDiameter_("diameter",centerMass_),
	fluidForce_("fluidForce",centerMass_),
	fluidTorque_("fluidTorque",centerMass_)
{

	auto domain = couplingMesh_.meshBox();

	if( 
		auto [domains, success] = processorComm_.collectAllToMaster(domain) ; 
		success)
	{
		
	}
	else
	{
		// - error 
		fatalErrorInFunction<<"could not corrlect meshBox into the master"<<endl;
		MPI::processor::abort(0);
	}
}


bool pFlow::coupling::couplingSystem::collectFluidForce()
{
	// realx3 scatteredComm is used 
	auto allForce = procDEMSystem_.particlesFluidForceAll();
	auto thisForce = makeSpan(fluidForce_);
	if(!realx3ScatteredComm_.collectSum(thisForce, allForce))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid force"<<endl;
		MPI::processor::abort(0);
	}

	return true;
}

bool pFlow::coupling::couplingSystem::collectFluidTorque()
{
	// realx3 scatteredComm is used 
	auto allTorque = procDEMSystem_.particlesFluidTorqueAll();
	auto thisTorque = makeSpan(fluidTorque_);
	if(!realx3ScatteredComm_.collectSum(thisTorque, allTorque))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid torque"<<endl;
		MPI::processor::abort(0);
	}

	return true;	
}

bool pFlow::coupling::couplingSystem::checkParticleDistribution()
{
	auto numParsInDomains = procDEMSystem_.numParInDomain();

	if( auto [thisNoPars, success] =  
		processorComm_.distributeMasterToAll(numParsInDomains); !success)
	{
		fatalErrorInFunction<<
		"failed to distribute particle numbers among processors"<<endl;
		MPI::processor::abort(0);
		return false;
	}
	else
	{
		if(!centerMass_.checkForNewSize(thisNoPars))
		{
			fatalErrorInFunction<<
			"cannot change the size of containers to new size "<< thisNoPars<<endl;
			MPI::processor::abort(0);
			return false;
		}
	}

	// first cunstructs index distribution
	auto parIndexInDomains = procDEMSystem_.parIndexInDomains();
	if(!realScatteredComm_.changeDataMaps(parIndexInDomains))
	{
		fatalErrorInFunction<<
		"error in creating index block for real type"<<endl;
		MPI::processor::abort(0);
		return false;
	}

	if(!realx3ScatteredComm_.changeDataMaps(parIndexInDomains))
	{
		fatalErrorInFunction<<
		"error in creating index block for realx3 type"<<endl;
		MPI::processor::abort(0);
		return false;
	}


	auto allDiam = procDEMSystem_.particledDiameterAll();
	auto thisDiam = makeSpan(particleDiameter_);

	if(!realScatteredComm_.distribute(allDiam, thisDiam))
	{
		fatalErrorInFunction<<
		"cannot distribute particle diameters to processors"<<endl;
		MPI::processor::abort(0);
		return false;
	}

	auto allPos = procDEMSystem_.particlesCenterMassAll();
	auto thisPos = makeSpan(centerMass_);
	if(!realx3ScatteredComm_.distribute(allPos, thisPos))
	{
		fatalErrorInFunction<<
		"cannot distribute particle positions to processors"<<endl;
		MPI::processor::abort(0);
		return false;
	}

	return true;

}