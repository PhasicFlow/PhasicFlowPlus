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

#ifndef __PCM_hpp__ 
#define __PCM_hpp__

#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{

/**
 * Particle In Cell (PCM) model for porosity calculation
 * 
 * This model only considers the particle center and if the particle center 
 * resides inside a cell, it is assumed that the whole volume of the particle
 * is located in that cell.
 * 
 */
class PCM
: 
	public porosity
{
public:

	/// Type info
	TypeInfo("PCM");

	/// Construc from dictionary 
	PCM(
		const unresolvedCouplingSystem& CS,
		const couplingMesh& 			cMesh,
		const Plus::realProcCMField& 	parDiam);

	/// Destructor
	virtual ~PCM() = default;

	/// Add this constructor to the list of virtual constructors
	add_vCtor
	(
		porosity,
		PCM,
		dictionary
	);

	bool internalFieldUpdate() override;


}; 

} // pFlow::coupling


#endif
