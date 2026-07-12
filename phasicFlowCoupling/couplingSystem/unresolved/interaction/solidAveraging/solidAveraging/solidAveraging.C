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

#include "solidAveraging.hpp"
#include "processorPlus.hpp"

pFlow::coupling::solidAveraging::solidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty,
    const word&                     name
)
:
    parAvFieldRef_(nullptr),
    uCS_(uCS),
    porosity_(prsty)
{
}


pFlow::uniquePtr<pFlow::coupling::solidAveraging> pFlow::coupling::solidAveraging::create
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty,
    const word&                     name
)
{
    
    if( wordvCtorSelector_.search(type))
	{
		Foam::Info<<"    Crearing solid averaging "
                  << Green_Text(type)
                  <<" for "<< Green_Text(name)<<"\n\n";
		return wordvCtorSelector_[type] (type, uCS, prsty, name);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor selector "<< type << " dose not exist"
				<<" for solid averaging method in " 
				<<"\nAvaiable ones are: \n"
				,
				wordvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}
