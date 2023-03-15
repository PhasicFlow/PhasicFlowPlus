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

#ifndef __subDivision9_hpp__ 
#define __subDivision9_hpp__

#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{


class subDivision9
: 
	public porosity
{
protected:

	int32 numInMesh_;

public:

	// type info
	TypeInfo("subDivision9");

	subDivision9(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	virtual ~subDivision9() = default;

	add_vCtor
	(
		porosity,
		subDivision9,
		dictionary
	);

	bool internalFieldUpdate() override;

	int32 numInMesh()const override
	{
		return numInMesh_;
	}

}; 

} // pFlow::coupling


#endif
