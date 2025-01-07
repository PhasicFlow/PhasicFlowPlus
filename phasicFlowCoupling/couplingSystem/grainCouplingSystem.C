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

#include "grainCouplingSystem.hpp"

bool pFlow::coupling::grainCouplingSystem::distributeParticleFields()
{
	

	auto allCG = pDEMSystem().particlesCourseGrainFactorMasterAllMaster();
	auto thisCG = makeSpan(courseGrainFactor_);
	if(!realMappedComm().distribute(allCG, thisCG))
	{
		fatalErrorInFunction<<
		"cannot distribute particle course grain factor among processors"<<endl;
		Plus::processor::abort(0);
		return false;
	}

	return couplingSystem::distributeParticleFields();
}

pFlow::coupling::grainCouplingSystem::grainCouplingSystem
(
	word demSystemName, 
	Foam::fvMesh& mesh,
	int argc, 
	char* argv[]
)
:
	couplingSystem(demSystemName, mesh, argc, argv),
	courseGrainFactor_("courseGrainFactor", centerMass()),
	porosityTimer_("porosity", &couplingTimers()),
	interactionTimer_("interaction", &couplingTimers())
{

	Foam::Info<<"\nCreating porosity model ...\n"<<Foam::endl;
	porosity_ = porosity::create(
		subDict("porosity"), 
		cMesh(),
		centerMass(),
		particleDiameter());

	Foam::Info<<"\nCreating grain drag model ...\n"<<Foam::endl;
	drag_ = grainDrag::create(
		subDict("drag"),
		porosity_()
		);
}


void pFlow::coupling::grainCouplingSystem::calculateFluidInteraction()
{
	interactionTimer_.start();
	if(drag_)
	{
		drag_->calculateGrainDragForce(
			particleVelocity(),
			particleDiameter(),
			courseGrainFactor_,
			fluidForce());
	}
	interactionTimer_.end();

	Foam::Info<<Blue_Text("Interaction time: ")<<interactionTimer_.lastTime()<<" s\n";
}

void pFlow::coupling::grainCouplingSystem::calculatePorosity()
{
	porosityTimer_.start();
	porosity_->calculatePorosity();
	porosity_->reportNumInMesh();
	porosityTimer_.end();

	Foam::Info<<Blue_Text("Porosity time: ")<<porosityTimer_.lastTime()<<" s\n";
}

const Foam::volScalarField& pFlow::coupling::grainCouplingSystem::alpha()const
{
	return porosity_->alpha();
}


const Foam::volScalarField& pFlow::coupling::grainCouplingSystem::Sp()const
{
	return drag_->Sp();
}


const Foam::volVectorField& pFlow::coupling::grainCouplingSystem::Su()const
{
	return drag_->Su();
}

