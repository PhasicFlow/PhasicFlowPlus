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


#include "PCM.hpp"
#include "self.hpp"


pFlow::coupling::PCM::PCM(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam)
:
	porosity(CS, cMesh, parDiam)
{

}

bool pFlow::coupling::PCM::internalFieldUpdate()
{
	
	self selfCellDist;
	
	auto solidVolTmp = calculateSolidVol(selfCellDist);

    Foam::fieldRef(*this) = Foam::max(
        1 - solidVolTmp/this->mesh().V(), 
        static_cast<Foam::scalar>(this->alphaMin()) );

	return true;
}
