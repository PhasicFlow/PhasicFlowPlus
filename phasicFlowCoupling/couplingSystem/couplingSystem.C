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
		word demSystemName, 
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
    MPI::procCommunication(),
	couplingMesh_(subDict("particleMapping"), mesh),
	procDEMSystem_(demSystemName, argc, argv),
	couplingTimers_("coupling", procDEMSystem_.getTimers()),
	cfdTimers_("CFD", procDEMSystem_.getTimers()),
	getDataTimer_("get data from DEM", &couplingTimers_),
	porosityTimer_("porosity", &couplingTimers_),
	interactionTimer_("interaction", &couplingTimers_),
	sendDataTimer_("send data to DEM", &couplingTimers_),
	centerMass_(),
	particleDiameter_("diameter",centerMass_),
	particleVelocity_("velocity", centerMass_),
	fluidForce_("fluidForce",centerMass_),
	fluidTorque_("fluidTorque",centerMass_)
{
	
	auto domain = couplingMesh_.meshBox();

	output<< " domain " << domain <<endl;

	if(! collectAllToAll(domain, meshBoxes_))
	{
		fatalErrorInFunction<<"could not corrlect meshBox into the master"<<endl;
		MPI::processor::abort(0);
	}

	// 
	Foam::Info<<"\nCreating porosity model ...\n"<<Foam::endl;

	porosity_ = porosity::create(
		subDict("porosity"), 
		couplingMesh_,
		centerMass_,
		particleDiameter_);

	Foam::Info<<"\nCreating drag model ...\n"<<Foam::endl;

	drag_ = drag::create(
		subDict("drag"),
		porosity_()
		);

}

bool pFlow::coupling::couplingSystem::getDataFromDEM(real t, real fluidDt)
{
	// this updates data on host side

	Foam::Info<<blueText("Obtaining data from DEM to master processor")<<Foam::endl;
	getDataTimer_.start();
	procDEMSystem_.getDataFromDEM();

	if( couplingMesh_.checkForDomainUpdate(t, fluidDt) )
	{
		Foam::Info<<blueText("Sub-domains have been updated at time ")<< 
		yellowText(t) <<Foam::endl;

		if(!procDEMSystem_.updateParticleDistribution(
			couplingMesh_.domainExpansionRatio(), 
			meshBoxes_))
		{
			fatalErrorInFunction;
			MPI::processor::abort(0);
			return false;
		}

		auto numParsInDomains = procDEMSystem_.numParInDomainMaster();

		if( auto [thisNoPars, success] =  
			distributeMasterToAll(numParsInDomains); !success)
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
				"cannot change the size of containers to new size "<< 
				thisNoPars<<endl;
				MPI::processor::abort(0);
				return false;
			}
			
		}
		
		// first cunstructs index distribution
		auto parIndexInDomains = procDEMSystem_.parIndexInDomainsMaster();

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
	}

	// update position and diameter in each processor
	distributeParticles();

	// update velocity in each processor
	distributeVelocity();

	getDataTimer_.end();
	return true;
}

bool pFlow::coupling::couplingSystem::sendDataToDEM()
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
	
	Foam::Info<<blueText("Sending fluid force from master processor to DEM")<<Foam::endl;

	if(!procDEMSystem_.sendFluidForceToDEM())
	{
		fatalErrorInFunction<< "could not perform sendFluidForceToDEM"<<endl;
		MPI::processor::abort(0);	
	}
}

void pFlow::coupling::couplingSystem::sendFluidTorqueToDEM()
{
	collectFluidTorque();
	
	Foam::Info<<blueText("Sending fluid torque from master processor to DEM")<<Foam::endl;

	if(!procDEMSystem_.sendFluidTorqueToDEM())
	{
		fatalErrorInFunction<< "could not perform sendFluidTorqueToDEM"<<endl;
		MPI::processor::abort(0);	
	}	
}

void pFlow::coupling::couplingSystem::calculateFluidInteraction()
{

	interactionTimer_.start();
	if(drag_)
	{
		drag_->calculateDragForce(
			particleVelocity_,
			particleDiameter_,
			fluidForce_);
	}
	interactionTimer_.end();

	output<<">>> Inteeraction time  "<<interactionTimer_.lastTime()<<endl;
	
}

void pFlow::coupling::couplingSystem::calculatePorosity()
{
	porosityTimer_.start();
	porosity_->calculatePorosity();
	porosity_->reportNumInMesh();
	porosityTimer_.end();
}



bool pFlow::coupling::couplingSystem::collectFluidForce()
{
	// realx3 scatteredComm is used 
	auto allForce = procDEMSystem_.particlesFluidForceAllMaster();
	for(size_t i=0; i<allForce.size(); i++)
		allForce[i] = zero3;

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
	auto allTorque = procDEMSystem_.particlesFluidTorqueAllMaster();
	for(size_t i=0; i<allTorque.size(); i++)
		allTorque[i] = zero3;
	auto thisTorque = makeSpan(fluidTorque_);
	if(!realx3ScatteredComm_.collectSum(thisTorque, allTorque))
	{
		fatalErrorInFunction<<
		"Faild to perform collective sum over processors for fluid torque"<<endl;
		MPI::processor::abort(0);
	}

	return true;	
}

bool pFlow::coupling::couplingSystem::distributeParticles()
{

	auto allDiam = procDEMSystem_.particledDiameterAllMaster();
	auto thisDiam = makeSpan(particleDiameter_);

	if(!realScatteredComm_.distribute(allDiam, thisDiam))
	{
		fatalErrorInFunction<<
		"cannot distribute particle diameters to processors"<<endl;
		MPI::processor::abort(0);
		return false;
	}

	auto allPos = procDEMSystem_.particlesCenterMassAllMaster();
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

bool pFlow::coupling::couplingSystem::distributeVelocity()
{
	auto allVel = procDEMSystem_.particlesVelocityAllMaster();
	auto thisVel = makeSpan(particleVelocity_);
	if(!realx3ScatteredComm_.distribute(allVel, thisVel))
	{
		fatalErrorInFunction<<
		"cannot distribute particle velocity among processors"<<endl;
		MPI::processor::abort(0);
		return false;
	}

	return true;
}