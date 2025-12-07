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

#include "fluidAveraging.hpp"
#include "unresolvedCouplingSystem.hpp"


pFlow::coupling::fluidAveraging::fluidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const word&                     name
)
:
    averagedField_
    (
        name,
        uCS.centerMass()
    ),
    uCS_(uCS)
{
}

pFlow::uniquePtr<pFlow::coupling::fluidAveraging> pFlow::coupling::fluidAveraging::create
(
    const word &type, 
    const unresolvedCouplingSystem &uCS, 
    const word &name
)
{
    if( wordvCtorSelector_.search(type))
	{
		REPORT(1) << "Creating fluid averaging of type "
				  << Green_Text(type) 
				  << " for " 
                  << Green_Text(name)
                  << END_REPORT;	
		return wordvCtorSelector_[type] (type, uCS, name);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< type << " dose not exist. \n\n"
				<<"\nAvaiable ones are: \n"
				,
				wordvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}
    return nullptr;
}
