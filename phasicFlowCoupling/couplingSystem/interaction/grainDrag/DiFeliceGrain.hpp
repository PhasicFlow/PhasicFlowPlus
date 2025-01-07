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

	Foam::scalar dimlessGrainDrag(Foam::scalar Re, Foam::scalar ep, Foam::scalar cgf)override
	{

		auto Rec = Foam::max(Re,residualRe_);
		Foam::scalar xi = 3.7 - 0.65*Foam::exp(-0.5*Foam::pow(1.5-Foam::log10(Rec),2));
		Foam::scalar Cd = Foam::pow(0.63+4.8/Foam::sqrt(Rec),2);

		Foam::scalar n = 0.014 + 0.955261 * Foam::exp( -4.10239227 * Foam::pow(ep, 0.002469095687*Rec)/Foam::pow(Rec, 0.41379176));
		n = Foam::max(Foam::min(n,1.0),0.0);

		return Cd/24 * Re * Foam::pow(ep, -xi)*Foam::pow(cgf,2-n); 
		
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
