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

#ifndef __DiFeliceGrain_hpp__ 
#define __DiFeliceGrain_hpp__

#include "grainDrag.hpp"

namespace pFlow::coupling
{
class DiFeliceGrain
:
	public grainDrag
{
protected:

	Foam::scalar dimlessGrainDrag(Foam::scalar Res, Foam::scalar ep, Foam::scalar cgf)override
	{

		auto Rec = Foam::max(Res,residualRe_);
		Foam::scalar xi = 3.7 - 0.65*Foam::exp(-0.5*Foam::pow(1.5-Foam::log10(Rec),2));
		Foam::scalar Cd = Foam::pow(0.63+4.8/Foam::sqrt(Rec),2);

		return Cd/24 * Res * Foam::pow(ep, -xi); 
		
	}

public:

	// type info
	TypeInfo("DiFeliceGrain");

	DiFeliceGrain(
		Foam::dictionary 		dict, 
		porosity& 				prsty);

	virtual ~DiFeliceGrain() = default;

	add_vCtor
	(
		grainDrag,
		DiFeliceGrain,
		dictionary
	);

	
}; 

} // pFlow::coupling


#endif
