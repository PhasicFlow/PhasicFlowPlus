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

#include "drag.hpp"
#include "unresolvedCouplingSystem.hpp"


void pFlow::coupling::drag::setSuSpToZero()
{
	// initialize all terms to zero
	forAll(Su_, celli)
	{
		Su_[celli] = Foam::vector(0,0,0);
		Sp_[celli] = 0.0;
	}
}

pFlow::coupling::drag::drag
(
	const unresolvedCouplingSystem& uCS, 
	const porosity& 				prsty
)
:
	Su_
	(
        Foam::IOobject
        (
            "Su",
            Foam::timeName(prsty.mesh().time()),
            prsty.mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::AUTO_WRITE
        ),
        prsty.mesh(),
        Foam::dimensionedVector
        (
            "Su",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
	),
    Sp_
    (
    	Foam::IOobject
	    (
	        "Sp",
	        Foam::timeName(prsty.mesh().time()),
	        prsty.mesh(),
	        Foam::IOobject::READ_IF_PRESENT,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	prsty.mesh(),
    	Foam::dimensionedScalar
    	(
    		"Sp", 
    		Foam::dimensionSet(1,-3,-1,0,0), 
    		0.0
    	)
    ),
	porosity_(prsty),
	p_
	(
		porosity_.mesh().lookupObject<Foam::volScalarField>("p")
	),
    isCompressible_( p_.dimensions() == Foam::dimPressure)
{}


Foam::tmp<Foam::volVectorField> 
pFlow::coupling::drag::pressureGradient(const Foam::volScalarField& rho)const
{
	if(isCompressible_)
		return Foam::fvc::grad(p_);
	else
		return Foam::fvc::grad(p_*rho);
}

const Foam::dictionary &pFlow::coupling::drag::dict() const
{
    return drag::getDict(porosity_.uCS());
}

const Foam::dictionary& pFlow::coupling::drag::getDict(const unresolvedCouplingSystem& uCS)
{
	return uCS.unresolvedDict().subDict("momentumInteraction").subDict("drag");
}

pFlow::uniquePtr<pFlow::coupling::drag> pFlow::coupling::drag::create
(
	const unresolvedCouplingSystem& uCS, 
	const porosity& 				prsty
)
{
	const auto& uDict = uCS.unresolvedDict();
	const auto& dragDict = drag::getDict(uCS);

	auto shapeName = uCS.shapeTypeName();
	auto dType = lookupDict<Foam::word>(dragDict, "model");
	
	Foam::word dragType;

	dragType = angleBracketsNames(shapeName+"Drag", dType);
	
	if( couplingSystemvCtorSelector_.search(dragType))
	{
		Foam::Info<<"    Crearing drag force "<<Green_Text(dragType)<<" ...\n\n";
		return couplingSystemvCtorSelector_[dragType] (uCS, prsty);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< dragType << " dose not exist"
				" for drag method in "<< uDict.name()
				<<"\nAvaiable ones are: \n"
				,
				couplingSystemvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}

