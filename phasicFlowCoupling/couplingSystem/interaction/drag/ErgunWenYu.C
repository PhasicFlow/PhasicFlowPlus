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


pFlow::coupling::ErgunWenYu::ErgunWenYu(
	Foam::dictionary 		dict, 
	porosity& 				prsty)
:
	drag(dict, prsty)
{

}

void pFlow::coupling::ErgunWenYu::calculateDragForce(
	const Foam::volVectorField& U,
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

	//const auto& mu = porosity_.mesh().lookupObject<Foam::volScalarField>("mu");
	//const auto& rho = porosity_.mesh().lookupObject<Foam::volScalarField>("rho");

	Foam::scalar mu = 1.8e-5;
	Foam::scalar rho = 1.2;

	const auto& parCells =  porosity_.particleCellIndex();


	for(size_t i=0; i<parCells.size(); i++)
	{
		auto cell = parCells[i];

		if(cell >= 0 )
		{
			auto Uf = U[cell];
			auto ef = porosity_[cell];

			Foam::vector up = {velocity[i].x(), velocity[i].y(),velocity[i].z()};
			Foam::scalar dp = diameter[i];
			
			Foam::vector ur = Uf-up;
			Foam::scalar Re = Foam::max(ef * rho * Foam::mag(ur) * dp /mu, residualRe_);

			Foam::scalar sp = 3 * Foam::constant::mathematical::pi * mu * ef * dp * dimlessDrag(Re, ef, dp);
			
			particleForce[i] = static_cast<real>(sp) * realx3(ur.x(), ur.y(), ur.z());
			Su_[cell] += -sp*up;
			Sp_[cell] += sp;

		}
	}

}

