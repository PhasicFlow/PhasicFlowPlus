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

#include "porosity.hpp"
#include "processor.hpp"

pFlow::coupling::porosity::porosity(
	Foam::dictionary		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	Foam::volScalarField
	(
		Foam::IOobject
	    (
	        "alpha",
	        cMesh.mesh().time().timeName(),
	        cMesh.mesh(),
	        Foam::IOobject::READ_IF_PRESENT,
	        Foam::IOobject::AUTO_WRITE
	    ),
    cMesh.mesh(),
    Foam::dimless
	),
	alphaMin_(dict.lookup<real>("alphaMin")),
	cMesh_(cMesh),
	centerMass_(centerMass),
	particleDiameter_(parDiam),
	parCellIndex_(
		"parCellIndex",
		static_cast<Foam::label>(-1),
		centerMass,
		true)
{

}


pFlow::uniquePtr<pFlow::coupling::porosity> 
	pFlow::coupling::porosity::create(
		Foam::dictionary		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam)
{
	auto method = dict.lookup<Foam::word>("method");
	if( dictionaryvCtorSelector_.search(method))
	{
		return dictionaryvCtorSelector_[method] (dict, cMesh, centerMass, parDiam);
	}
	else
	{
		if(MPI::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< method << " dose not exist"
				" for porosity method in "<< dict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		MPI::processor::abort(0);
	}

	return nullptr;
}