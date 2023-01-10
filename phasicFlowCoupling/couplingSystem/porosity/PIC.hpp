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

#ifndef __PIC_hpp__ 
#define __PIC_hpp__

#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{


class PIC
: 
	public porosity
{
protected:

	

public:

	// type info
	TypeInfo("PIC");

	PIC(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	virtual ~PIC() = default;

	add_vCtor
	(
		porosity,
		PIC,
		dictionary
	);

	bool calculatePorosity() override;

}; 

} // pFlow::coupling


#endif
