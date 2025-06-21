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
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[],
		bool requireRVel)
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
	particleMapping_(subDict("particleMapping")),
	couplingMesh_
	(
		subDict("particleMapping"), 
		mesh,
		particleMapping_.centerMass()
	),
	procDEMSystem_
	(
		shapeTypeName+"DEMSystem", 
		argc, 
		argv, 
		requireRVel
	),
	couplingTimers_
	(
		"coupling", 
		procDEMSystem_.getTimers()
	),
	cfdTimers_
	(
		"CFD", 
		procDEMSystem_.getTimers()
	),
	getDataTimer_
	(
		"get data from DEM", 
		&couplingTimers_
	),
	sendDataTimer_
	(
		"send data to DEM", 
		&couplingTimers_
	),
	particleID_
	(
		"particleID",
		particleMapping_.centerMass()
	),
	particleDiameter_
	(
		"diameter",
		particleMapping_.centerMass()
	),
	particleVelocity_
	(
		"velocity", 
		particleMapping_.centerMass()
	),
	particleRVelocity_(
		"rVelocity", 
		particleMapping_.centerMass()
	),
	fluidForce_(
		"fluidForce",
		particleMapping_.centerMass()
	),
	fluidTorque_(
		"fluidTorque",
		particleMapping_.centerMass()
	),
	requireRVel_(requireRVel)
{}

bool pFlow::coupling::couplingSystem::getDataFromDEM(real t, real fluidDt)
{
	// this updates data on host side

	Foam::Info<<Blue_Text("Obtaining data from DEM to master processor")<<Foam::endl;
	getDataTimer_.start();
	procDEMSystem_.getDataFromDEM();

	if( !particleMapping_.update(t, fluidDt, procDEMSystem_, couplingMesh_) ) return false;

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
	
	if(!particleMapping_.realx3ScatteredComm().collectSum(thisForce, allForce))
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
	if(!particleMapping_.realx3ScatteredComm().collectSum(thisTorque, allTorque))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid torque"<<endl;
		Plus::processor::abort(0);
	}

	return true;	
}


bool pFlow::coupling::couplingSystem::distributeParticleFields()
{
	
    auto allDiam = procDEMSystem_.particlesDiameterAllMaster();
    auto thisDiam = makeSpan(particleDiameter_);

    if(!particleMapping_.realScatteredComm().distribute(allDiam, thisDiam))
    {
        fatalErrorInFunction<<
        "cannot distribute particle diameters to processors"<<endl;
        Plus::processor::abort(0);
        return false;
    }

    auto allVel = procDEMSystem_.particlesVelocityAllMaster();
    auto thisVel = makeSpan(particleVelocity_);
    if(!particleMapping_.realx3ScatteredComm().distribute(allVel, thisVel))
    {
        fatalErrorInFunction<<
        "cannot distribute particle velocity among processors"<<endl;
        Plus::processor::abort(0);
        return false;
    }

    if( requireRVel_)
    {
        auto allRVel = procDEMSystem_.particlesRVelocityAllMaster();
        auto thisRVel = makeSpan(particleRVelocity_);
        if(!particleMapping_.realx3ScatteredComm().distribute(allRVel, thisRVel))
        {
            fatalErrorInFunction<<
            "cannot distribute particle rotational velocity among processors"<<endl;
            Plus::processor::abort(0);
            return false;
        }
    }

    /*auto allID = procDEMSystem_.particleIdAllMaster();
    auto thisID = makeSpan(particleID_);

    if(!particleMapping_.uint32ScatteredComm().distribute(allID, thisID))
    {
        fatalErrorInFunction<<
        "cannot distribute particle IDs to processors"<<endl;
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
