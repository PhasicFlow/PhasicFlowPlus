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
// from OpenFOAM
#include "mathematicalConstants.H"
#include "volFields.H"

#include "ErgunWenYu.hpp"
#include "streams.hpp"
#include "processor.hpp"


pFlow::coupling::ErgunWenYu::ErgunWenYu(
	Foam::dictionary 		dict, 
	porosity& 				prsty)
:
	drag(dict, prsty)
{

}

void pFlow::coupling::ErgunWenYu::calculateDragForce(
	const MPI::realx3ProcCMField& velocity,
	const MPI::realProcCMField& diameter,
	MPI::realx3ProcCMField& particleForce)
{

	// initialize all terms to zero
	Su_ = Foam::dimensionedVector(
	    	"Su",
	    	Foam::dimensionSet(1,-2,-2,0,0),
	    	Foam::vector(0,0,0));

	Sp_ = Foam::dimensionedScalar(
    		"Sp", 
    		Foam::dimensionSet(1,-3,-1,0,0), 
    		0.0);
	particleForce = realx3(0,0,0);

	// gets pressure gradient 
	auto pGrad = pressureGradient();
	auto& pGradRef = pGrad();

	const auto& parCells =  porosity_.particleCellIndex();
	if(parCells.size() == 0 ) return;

	const auto& nu = porosity_.mesh().lookupObject<Foam::volScalarField>("nu");
	const auto& rho = porosity_.mesh().lookupObject<Foam::volScalarField>("rho");

	for(size_t i=0; i<parCells.size(); i++)
	{
		auto cell = parCells[i];

		if(cell >= 0 )
		{
			auto rhoi = rho[cell];
			auto mui = nu[cell]* rhoi;
			auto Uf = U_[cell];
			auto ef = porosity_[cell];
			auto dp = diameter[i];

			Foam::vector up = {velocity[i].x(), velocity[i].y(),velocity[i].z()};

			auto vp = Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);
			
			Foam::vector ur = Uf-up;
			Foam::scalar Re = Foam::max(ef * rhoi * Foam::mag(ur) * dp /mui, residualRe_);

			Foam::scalar sp = 3 * Foam::constant::mathematical::pi * mui * ef * dp * dimlessDrag(Re, ef, dp);
			
			Foam::vector pf = {0.00000015, 0, 0}; // static_cast<real>(sp)*ur; //- vp*pGradRef[cell];
			particleForce[i] = realx3(pf.x(), pf.y(), pf.z());
			Su_[cell] += -sp*up;
			Sp_[cell] += sp;
		}
	}

}