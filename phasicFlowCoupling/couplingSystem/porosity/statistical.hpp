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

#ifndef __statistical_hpp__ 
#define __statistical_hpp__

#include "virtualConstructor.hpp"
// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{

/**
 * Particle In Cell (statistical) model for porosity calculation
 * 
 * This model only considers the particle center and if the particle center 
 * resides inside a cell, it is assumed that the whole volume of the particle
 * is located in that cell.
 * 
 */
class statistical
: 
	public porosity
{
protected:

	int32 numInMesh_ 		= 0;


	// search/
	Foam::scalar 				b_;

	std::vector<std::vector<Foam::label>>  neighborList_;		

public:

	/// Type info
	TypeInfo("statistical");

	/// Construc from dictionary 
	statistical(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	/// Destructor
	virtual ~statistical() = default;

	/// Add this constructor to the list of virtual constructors
	add_vCtor
	(
		porosity,
		statistical,
		dictionary
	);

	bool internalFieldUpdate() override;


}; 

} // pFlow::coupling


#endif
