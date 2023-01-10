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

bool pFlow::coupling::couplingSystem::checkForDomainUpdate(real t, real fluidDt)
{
	real sT = procDEMSystem_.startTime();
	if( pFlow::equal(t, sT) ) 
	{
		lastTimeUpdated_ = t;
		return true;
	}

	if( abs(t-(lastTimeUpdated_+subDomainUpdateInterval_)) < 0.98*fluidDt)
	{
		lastTimeUpdated_ = t;
		return true;
	}

	return false;

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
            mesh.time().caseConstant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    subDomainExpansionFraction_(lookup<pFlow::real>("subDomainExpansionFraction")),
    subDomainUpdateInterval_(lookup<pFlow::real>("subDomainUpdateInterval")),
	couplingMesh_(mesh),
	processorComm_(),
	procDEMSystem_(demSystemName, argc, argv),
	centerMass_(),
	particleDiameter_("diameter",centerMass_),
	particleVelocity_("velocity", centerMass_),
	fluidForce_("fluidForce",centerMass_),
	fluidTorque_("fluidTorque",centerMass_)
{
	auto domain = couplingMesh_.meshBox();

	if(! processorComm_.collectAllToAll(domain, meshBoxes_))
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
	procDEMSystem_.getDataFromDEM();

	
	if( checkForDomainUpdate(t, fluidDt) )
	{
		Foam::Info<<"Sub-domains have been updated at time "<< t<<Foam::endl;

		if(!procDEMSystem_.updateParticleDistribution(
			 subDomainExpansionFraction_, 
			meshBoxes_))
		{
			fatalErrorInFunction;
			MPI::processor::abort(0);
			return false;
		}

		auto numParsInDomains = procDEMSystem_.numParInDomainMaster();

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

	return true;
}

bool pFlow::coupling::couplingSystem::calculatePorosity()
{
	return porosity_->calculatePorosity();
}

bool pFlow::coupling::couplingSystem::collectFluidForce()
{
	// realx3 scatteredComm is used 
	auto allForce = procDEMSystem_.particlesFluidForceAllMaster();
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