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
#include "volFields.H"
#include "fvc.H"

#include "grainDrag.hpp"
#include "processorPlus.hpp"


pFlow::coupling::grainDrag::grainDrag(
	Foam::dictionary 		dict, 
	porosity& 				prsty)
:
	drag(dict, prsty)
{

}


void pFlow::coupling::grainDrag::calculateGrainDragForce(
	const Plus::realx3ProcCMField& velocity,
	const Plus::realProcCMField& diameter,
	const Plus::realProcCMField& courseGrainFactor,
	Plus::realx3ProcCMField& particleForce)
{

	setSuSpToZero();
	particleForce = realx3(0,0,0);


	const auto& parCells =  porosity_.particleCellIndex();

	const auto& nu = porosity_.mesh().lookupObject<Foam::volScalarField>("nu");
	const auto& rho = porosity_.mesh().lookupObject<Foam::volScalarField>("rho");

	// gets pressure gradient 
	auto pGrad = pressureGradient(rho);
	auto& pGradRef = pGrad();
	const auto& Vcell = Su_.mesh().V();

	size_t numPar = parCells.size();

#pragma omp parallel for
	for(size_t i=0; i<numPar; i++)
	{
		auto cell = parCells[i];

		if(cell >= 0 )
		{
			auto rhoi = rho[cell];
			auto mui = nu[cell]* rhoi;
			auto Uf = U_[cell];
			auto ef = porosity_[cell];
			auto cgf = courseGrainFactor[i];
			auto dp = diameter[i];
			auto dps = dp/cgf; // dps = orginal particles diameter

			Foam::vector up = {velocity[i].x(), velocity[i].y(),velocity[i].z()};

			auto vp = Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);
			
			Foam::vector ur = Uf-up;
			Foam::scalar Res = Foam::max(ef * rhoi * Foam::mag(ur) * dps /mui, residualRe_) ;

			Foam::scalar sp = 3*Foam::pow(cgf,3)*Foam::constant::mathematical::pi * mui * ef * dps * dimlessGrainDrag(Res, ef, cgf);
			
			Foam::vector pf = static_cast<real>(sp)*ur - vp*pGradRef[cell];
			
			
			particleForce[i] = realx3(pf.x(), pf.y(), pf.z());
			
			Foam::vector suu = (-sp*up)/Vcell[cell];
			#pragma omp atomic
			Su_[cell].x() += suu.x(); 

			#pragma omp atomic
			Su_[cell].y() += suu.y(); 

			#pragma omp atomic
			Su_[cell].z() += suu.z(); 
			
			#pragma omp atomic
			Sp_[cell] += sp/Vcell[cell];
		}
	}

}

pFlow::uniquePtr<pFlow::coupling::grainDrag> 
pFlow::coupling::grainDrag::create(
	Foam::dictionary 		dict, 
	porosity& 				prsty)
{
	auto type = dict.lookup<Foam::word>("type");
	if( dictionaryvCtorSelector_.search(type))
	{
		return dictionaryvCtorSelector_[type] (dict, prsty);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< type << " dose not exist"
				" for grain drag method in "<< dict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}