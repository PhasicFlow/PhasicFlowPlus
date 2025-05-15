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


bool pFlow::coupling::couplingSystem::updateMeshBoxes()
{

	auto domain = couplingMesh_.meshBox();

	if(! collectAllToAll(domain, meshBoxes_))
	{
		fatalErrorInFunction<<"could not corrlect meshBox into the master"<<endl;
		return false;
	}

	return true;
}

pFlow::coupling::couplingSystem::couplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[])
:
	Foam::IOdictionary
    (
        IOobject
        (
            "couplingProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Plus::procCommunication(),
	couplingMesh_(subDict("particleMapping"), mesh),
	procDEMSystem_(shapeTypeName+"DEMSystem", argc, argv),
	couplingTimers_("coupling", procDEMSystem_.getTimers()),
	cfdTimers_("CFD", procDEMSystem_.getTimers()),
	getDataTimer_("get data from DEM", &couplingTimers_),
	sendDataTimer_("send data to DEM", &couplingTimers_),
	centerMass_(),
	particleID_("particleID",centerMass_),
	particleDiameter_("diameter",centerMass_),
	particleVelocity_("velocity", centerMass_),
	particleRVelocity_("rVelocity", centerMass_),
	fluidForce_("fluidForce",centerMass_),
	fluidTorque_("fluidTorque",centerMass_)
{
	
	auto domain = couplingMesh_.meshBox();

	if(! collectAllToAll(domain, meshBoxes_))
	{
		fatalErrorInFunction<<"could not corrlect meshBox into the master"<<endl;
		Plus::processor::abort(0);
	}

	pFlow::mOutput<<"meshBoxes:\n"<<meshBoxes_<<pFlow::endl;

}

bool pFlow::coupling::couplingSystem::getDataFromDEM(real t, real fluidDt)
{
	// this updates data on host side

	Foam::Info<<Blue_Text("Obtaining data from DEM to master processor")<<Foam::endl;
	getDataTimer_.start();
	procDEMSystem_.getDataFromDEM();

	if( couplingMesh_.checkForDomainUpdate(t, fluidDt) )
	{
		Foam::Info<<Blue_Text("Sub-domains have been updated at time ")<< 
		Yellow_Text(t) <<Foam::endl;

		if(!procDEMSystem_.updateParticleDistribution(
			couplingMesh_.domainExpansionRatio(), 
			meshBoxes_))
		{
			fatalErrorInFunction;
			Plus::processor::abort(0);
			return false;
		}

		auto numParsInDomains = procDEMSystem_.numParInDomainMaster();

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
		auto parIndexInDomains = procDEMSystem_.parIndexInDomainsMaster();

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
	}

	// update position and diameter in each processor
	distributeParticles();

	// update velocity in each processor
	distributeParticleFields();

	getDataTimer_.end();
	return true;
}

bool pFlow::coupling::couplingSystem::sendDataToDEM(real, real)
{
	sendDataTimer_.start();
		sendFluidForceToDEM();
		sendFluidTorqueToDEM();
	sendDataTimer_.end();
	return true;
}

void pFlow::coupling::couplingSystem::sendFluidForceToDEM()
{
	collectFluidForce();
	
	Foam::Info<<Blue_Text("Sending fluid force from master processor to DEM")<<Foam::endl;

	if(!procDEMSystem_.sendFluidForceToDEM())
	{
		fatalErrorInFunction<< "could not perform sendFluidForceToDEM"<<endl;
		Plus::processor::abort(0);	
	}
}

void pFlow::coupling::couplingSystem::sendFluidTorqueToDEM()
{
	collectFluidTorque();
	
	Foam::Info<<Blue_Text("Sending fluid torque from master processor to DEM")<<Foam::endl;

	if(!procDEMSystem_.sendFluidTorqueToDEM())
	{
		fatalErrorInFunction<< "could not perform sendFluidTorqueToDEM"<<endl;
		Plus::processor::abort(0);	
	}	
}



bool pFlow::coupling::couplingSystem::collectFluidForce()
{
	// realx3 scatteredComm is used 
	auto allForce = procDEMSystem_.particlesFluidForceAllMaster();
	
	for(uint32 i=0; i<allForce.size(); i++)
		allForce[i] = zero3;

	auto thisForce = makeSpan(fluidForce_);
	
	if(!realx3ScatteredComm_.collectSum(thisForce, allForce))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid force"<<endl;
		Plus::processor::abort(0);
	}

	return true;
}

bool pFlow::coupling::couplingSystem::collectFluidTorque()
{
	// realx3 scatteredComm is used 
	auto allTorque = procDEMSystem_.particlesFluidTorqueAllMaster();
	for(uint32 i=0; i<allTorque.size(); i++)
		allTorque[i] = zero3;
	auto thisTorque = makeSpan(fluidTorque_);
	if(!realx3ScatteredComm_.collectSum(thisTorque, allTorque))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid torque"<<endl;
		Plus::processor::abort(0);
	}

	return true;	
}

bool pFlow::coupling::couplingSystem::distributeParticles()
{
	auto allDiam = procDEMSystem_.particlesDiameterAllMaster();
	auto thisDiam = makeSpan(particleDiameter_);

	if(!realScatteredComm_.distribute(allDiam, thisDiam))
	{
		fatalErrorInFunction<<
		"cannot distribute particle diameters to processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	auto allPos = procDEMSystem_.particlesCenterMassAllMaster();
	auto thisPos = makeSpan(centerMass_);
	if(!realx3ScatteredComm_.distribute(allPos, thisPos))
	{
		fatalErrorInFunction<<
		"cannot distribute particle positions to processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	/*auto allID = procDEMSystem_.particlesIDAllMaster();
	auto thisID = makeSpan(particleID_);

	if(!realScatteredComm_.distribute(allID, thisID))
	{
		fatalErrorInFunction<<
		"cannot distribute particle IDs to processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}*/

	return true;

}

bool pFlow::coupling::couplingSystem::distributeParticleFields()
{
	auto allVel = procDEMSystem_.particlesVelocityAllMaster();
	auto thisVel = makeSpan(particleVelocity_);
	if(!realx3ScatteredComm_.distribute(allVel, thisVel))
	{
		fatalErrorInFunction<<
		"cannot distribute particle velocity among processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	/*auto allRVel = procDEMSystem_.particlesRVelocityAllMaster();
	auto thisRVel = makeSpan(particleRVelocity_);
	if(!realx3ScatteredComm_.distribute(allRVel, thisRVel))
	{
		fatalErrorInFunction<<
		"cannot distribute particle rotational velocity among processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}*/

	return true;
}

bool pFlow::coupling::couplingSystem::iterate(real upToTime, bool writeTime, const word& timeName)
{
	Foam::Info<<Blue_Text("Iterating DEM upto time ") << Yellow_Text(upToTime)<<Foam::endl;
	return procDEMSystem_.iterate(upToTime, writeTime, timeName);
}