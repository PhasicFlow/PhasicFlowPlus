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

#include "drag.hpp"
#include "processor.hpp"

pFlow::coupling::drag::drag(
	Foam::dictionary 		dict, 
	porosity& 					prsty)
:
	residualRe_(dict.lookup<real>("residualRe")),
	porosity_(prsty),
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
    Foam::dimensionSet(1,-2,-2,0,0)),
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
    	Foam::dimensionSet(1,-3,-1,0,0))
{

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
		if(MPI::processor::isMaster())
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
		MPI::processor::abort(0);
	}

	return nullptr;
}