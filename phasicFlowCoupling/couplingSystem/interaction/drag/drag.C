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

#include "drag.hpp"
#include "processorPlus.hpp"

void pFlow::coupling::drag::setSuSpToZero()
{
	// initialize all terms to zero
	Su_ = Foam::dimensionedVector(
    	Su_.name(),
    	Su_.dimensions(),
    	Foam::vector(0,0,0));

	Sp_ = Foam::dimensionedScalar(
		Sp_.name(), 
		Sp_.dimensions(), 
		0.0);
}

pFlow::coupling::drag::drag(
	Foam::dictionary 		dict, 
	porosity& 				prsty)
:
	residualRe_(dict.lookup<real>("residualRe")),
	porosity_(prsty),
	p_(
		porosity_.mesh().lookupObject<Foam::volScalarField>("p")
	),
	U_(
		porosity_.mesh().lookupObject<Foam::volVectorField>("U")
		),
	Su_(
	    Foam::IOobject
	    (
	        "Su",
	        porosity_.mesh().time().timeName(),
	        porosity_.mesh(),
	        Foam::IOobject::READ_IF_PRESENT,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	porosity_.mesh(),
    	Foam::dimensionedVector(
	    	"Su",
	    	Foam::dimensionSet(1,-2,-2,0,0),
	    	Foam::vector(0,0,0))
		),
    Sp_(
    	Foam::IOobject
	    (
	        "Sp",
	        porosity_.mesh().time().timeName(),
	        porosity_.mesh(),
	        Foam::IOobject::READ_IF_PRESENT,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	porosity_.mesh(),
    	Foam::dimensionedScalar(
    		"Sp", 
    		Foam::dimensionSet(1,-3,-1,0,0), 
    		0.0)
    	),
    isCompressible_(
    	p_.dimensions() == Foam::dimPressure
    	)
{

}

Foam::tmp<Foam::volVectorField> 
pFlow::coupling::drag::pressureGradient(const Foam::volScalarField& rho)const
{
	if(isCompressible_)
		return Foam::fvc::grad(p_);
	else
		return Foam::fvc::grad(p_)*rho;
}


void pFlow::coupling::drag::calculateDragForce(
	const Plus::realx3ProcCMField& velocity,
	const Plus::realProcCMField& diameter,
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
			auto dp = diameter[i];

			Foam::vector up = {velocity[i].x(), velocity[i].y(),velocity[i].z()};

			auto vp = Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);
			
			Foam::vector ur = Uf-up;
			Foam::scalar Re = Foam::max(ef * rhoi * Foam::mag(ur) * dp /mui, residualRe_);

			Foam::scalar sp = 3 * Foam::constant::mathematical::pi * mui * ef * dp * dimlessDrag(Re, ef);
			
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

pFlow::uniquePtr<pFlow::coupling::drag> 
	pFlow::coupling::drag::create(
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
				" for porosity method in "<< dict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}