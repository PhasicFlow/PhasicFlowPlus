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

#ifndef __grainDrag_hpp__ 
#define __grainDrag_hpp__


// from phasicFlow-coupling
#include "drag.hpp"

namespace pFlow::coupling
{
class grainDrag
:
	public drag
{
protected:

	
	virtual 
	Foam::scalar dimlessGrainDrag(Foam::scalar Re, Foam::scalar ep, Foam::scalar cgf)=0;

	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)override 
	{
		notImplementedFunction;
		return 0;
	}

	using drag::calculateDragForce;
public:

	// type info
	TypeInfo("grainDrag");

	grainDrag(
		Foam::dictionary 		dict, 
		porosity& 				prsty);

	virtual ~grainDrag() = default;

	create_vCtor
	(
		grainDrag,
		dictionary,
		(
			Foam::dictionary 		dict, 
			porosity& 				prsty
		),
		(dict, prsty)
	);

	void calculateGrainDragForce(
		const Plus::realx3ProcCMField& velocity,
		const Plus::realProcCMField& diameter,
		const Plus::realProcCMField& courseGrainFactor,
		Plus::realx3ProcCMField& particleForce);



	static
	uniquePtr<grainDrag> create(
		Foam::dictionary 		dict, 
		porosity& 				prsty);
		
	
}; 

} // pFlow::coupling


#endif
