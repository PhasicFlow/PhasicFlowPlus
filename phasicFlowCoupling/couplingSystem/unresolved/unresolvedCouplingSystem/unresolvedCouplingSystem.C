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

#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::unresolvedCouplingSystem::unresolvedCouplingSystem
(
    word demSystemName,
    word couplingSystemType, 
    Foam::fvMesh& mesh,
    int argc, 
    char* argv[]
)
:
    couplingSystem(demSystemName, mesh, argc, argv, true)
{
    distribution_ = distributionBase::create
    (
        unresolvedDict(),
        cMesh(),
        parMapping().centerMass()
    );
}


const Foam::dictionary& pFlow::coupling::unresolvedCouplingSystem::unresolvedDict()const
{
    return this->subDict("unresolved");
}


pFlow::uniquePtr<pFlow::coupling::unresolvedCouplingSystem> 
    pFlow::coupling::unresolvedCouplingSystem::create
    (
        word shapeTypeName,
        word couplingSystemType,
        Foam::fvMesh& mesh,
        int argc, 
        char* argv[]
    )
{

    word couplingType = angleBracketsNames(
    	shapeTypeName +"UnresolvedCouplingSystem",
		couplingSystemType);

    if( wordvCtorSelector_.search(couplingType))
    {
        REPORT(0) << "Creating unresolved coupling system "
                    << Green_Text(couplingType) 
                    << " ..." << END_REPORT;	
        return wordvCtorSelector_[couplingType] (shapeTypeName, couplingSystemType, mesh, argc, argv);
    }
    else
    {
        if(Plus::processor::isMaster())
        {
            printKeys
            ( 
                fatalErrorInFunction << "Ctor Selector "<< couplingType << " dose not exist. \n\n"
                <<"\nAvaiable ones are: \n"
                ,
                wordvCtorSelector_
            )<<endl;
        }
        Plus::processor::abort(0);
    }

	return nullptr;

}

Foam::tmp<Foam::volVectorField> 
    pFlow::coupling::unresolvedCouplingSystem::Us() const
{
    // Default fallback implementation for the velocity of injected mass source.
    // Returns a continuous zero vector field mapped over the fluid domain.
    // Safely supports physical setups where mass injection holds zero initial 
    // momentum or when spatial distribution is bypassed.
    return Foam::tmp<Foam::volVectorField>::New(
        Foam::IOobject(
            "UsZero",
            this->cMesh().mesh().time().timeName(),
            this->cMesh().mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE,
            false),
        this->cMesh().mesh(),
        Foam::dimensionedVector(
            "zero", 
            Foam::dimVelocity, 
            Foam::vector::zero));
}