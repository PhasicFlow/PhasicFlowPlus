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


template<typename DistributorType>
pFlow::coupling::sphereUnresolvedCouplingSystem<DistributorType>::sphereUnresolvedCouplingSystem
(
	word demSystemName, 
	Foam::fvMesh& mesh,
	int argc, 
	char* argv[]
)
:
	UnresolvedCouplingSystem<DistributorType>(demSystemName, mesh, argc, argv),
	porosityTimer_("porosity", &this->couplingTimers()),
	interactionTimer_("interaction", &this->couplingTimers())
{

	Foam::Info<<"\n  Creating porosity model ..."<<Foam::endl;
	porosity_ = porosity::create(
		*this, 
		this->cMesh(),
		this->particleDiameter());

	Foam::Info<<"\n  Creating drag model ..."<<Foam::endl;
	drag_ = drag<DistributorType>::create(
		*this,
		porosity_()
		);

	requiresDistribution_ = 
		porosity_().requireCellDistribution()|| drag_().requireCellDistribution();
}

template<typename DistributorType>
void pFlow::coupling::sphereUnresolvedCouplingSystem<DistributorType>::calculatePorosity()
{
	porosityTimer_.start();

	// update coupling mesh and map particles 
	this->cMesh().update();
	
	// update weights for distribution (if coupling requires it)
	if(requiresDistribution_)
		this->cellDistribution().updateWeights(this->cMesh().parCellIndex(), this->particleDiameter());

	// calculate porosity 
	porosity_->calculatePorosity();

	porosityTimer_.end();

	Foam::Info<<Blue_Text("Porosity time: ")<< 
				Yellow_Text(porosityTimer_.lastTime())<<
				Yellow_Text(" s")<<Foam::endl;
}

template<typename DistributorType>
void pFlow::coupling::sphereUnresolvedCouplingSystem<DistributorType>::calculateFluidInteraction()
{
	const auto& U = this->cMesh().mesh().template lookupObject<Foam::volVectorField>("U");
	interactionTimer_.start();
	if(drag_)
	{
		drag_->calculateDragForce(
			U,
			this->particleVelocity(),
			this->particleDiameter(),
			this->fluidForce());
	}
	interactionTimer_.end();

	Foam::Info<<Blue_Text("Interaction time: ")<< 
				Yellow_Text(interactionTimer_.lastTime())<<
				Yellow_Text(" s")<<Foam::endl;
}
