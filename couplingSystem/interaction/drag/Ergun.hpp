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

#ifndef __Ergun_hpp__ 
#define __Ergun_hpp__

#include "drag.hpp"

namespace pFlow::coupling
{
class Ergun
:
	public drag
{
protected:

public:

	// type info
	TypeInfo("Ergun");

	Ergun(
		Foam::dictionary 		dict, 
		porosity& 				prsty);

	virtual ~Ergun() = default;

	add_vCtor
	(
		drag,
		Ergun,
		dictionary
	);

	bool calculateDragForce(
		const Foam::volVectorField& U,
		const MPI::realx3ProcCMField& velocity,
		const MPI::realProcCMField& diameter,
		MPI::realx3ProcCMField& particleForce)override;
		
	
}; 

} // pFlow::coupling


#endif
