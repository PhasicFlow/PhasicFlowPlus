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


/**
 * @class noneDrag
 * @brief This is a drag force model to disable drag calculation in the coupling system.
 * 
 * 
 * @see dimlessDrag()
 */
#ifndef __none_hpp__ 
#define __none_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"


#include "typeInfo.hpp"

namespace pFlow::coupling
{

class noneDrag
{

public:

	// type info
	TypeInfoNV("none");

	noneDrag(const Foam::dictionary& dict);

	inline
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep) const
	{
		return 0.0;
	}
	
	inline 
	Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
	{
		return dimlessDrag(Re, ep);
	}
	
}; 

} // pFlow::coupling


#endif
