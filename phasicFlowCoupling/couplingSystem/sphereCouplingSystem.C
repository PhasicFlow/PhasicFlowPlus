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

#include "sphereCouplingSystem.hpp"

pFlow::coupling::sphereCouplingSystem::sphereCouplingSystem
(
	word demSystemName, 
	Foam::fvMesh& mesh,
	int argc, 
	char* argv[]
)
:
	couplingSystem(demSystemName, mesh, argc, argv),
	porosityTimer_("porosity", &couplingTimers()),
	interactionTimer_("interaction", &couplingTimers())
{

	Foam::Info<<"\nCreating porosity model ...\n"<<Foam::endl;
	porosity_ = porosity::create(
		subDict("porosity"), 
		cMesh(),
		centerMass(),
		particleDiameter());

	Foam::Info<<"\nCreating drag model ...\n"<<Foam::endl;
	drag_ = drag::create(
		subDict("drag"),
		porosity_()
		);
}


void pFlow::coupling::sphereCouplingSystem::calculateFluidInteraction()
{
	interactionTimer_.start();
	if(drag_)
	{
		drag_->calculateDragForce(
			particleVelocity(),
			particleDiameter(),
			fluidForce());
	}
	interactionTimer_.end();

	Foam::Info<<Blue_Text("Interaction time: ")<<interactionTimer_.lastTime()<<" s\n";
}

void pFlow::coupling::sphereCouplingSystem::calculatePorosity()
{
	porosityTimer_.start();
	porosity_->calculatePorosity();
	porosity_->reportNumInMesh();
	porosityTimer_.end();

	Foam::Info<<Blue_Text("Porosity time: ")<<porosityTimer_.lastTime()<<" s\n";
}

const Foam::volScalarField& pFlow::coupling::sphereCouplingSystem::alpha()const
{
	return porosity_->alpha();
}


const Foam::volScalarField& pFlow::coupling::sphereCouplingSystem::Sp()const
{
	return drag_->Sp();
}


const Foam::volVectorField& pFlow::coupling::sphereCouplingSystem::Su()const
{
	return drag_->Su();
}

