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


// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{

/**
 * sub-division model for calculating fluid porosity
 * 
 * This model divides a sphere into 9 equal parts and locate each of these 
 * parts in the cells. The volume of each part located in each cell is 
 * considered as the solid volume in that cell. 
 * 
 */
class subDivision9
: 
	public porosity
{
public:

	/// Type info
	TypeInfo("subDivision9");

	/// Construct from dictionary
	subDivision9(
		const unresolvedCouplingSystem& CS,
		const couplingMesh& 			cMesh,
		const Plus::realProcCMField& 	parDiam);

	/// Destructor
	virtual ~subDivision9() = default;

	/// Add this constructor to the list of virtual constructors
	add_vCtor
	(
		porosity,
		subDivision9,
		dictionary
	);

	bool internalFieldUpdate() override;


}; 

} // pFlow::coupling


#endif
