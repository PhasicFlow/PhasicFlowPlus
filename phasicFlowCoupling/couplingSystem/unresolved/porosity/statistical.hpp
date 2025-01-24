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
 * statistical model for porosity calculation
 * 
 * This model distributes particle volume to the target cell and neighboring 
 * cell based on the a Gaussian distribution and the distance between 
 * particle center and cell center. 
 * 
 */
class statistical
: 
	public porosity
{
protected:

	// radius of circule for cell neighbor list 
	Foam::scalar 				neighborLength_;

	// list of nieghbor cells for each cell 
	std::vector<std::vector<Foam::label>>  neighborList_;
	
	// boundary cells 
	std::vector<std::pair<Foam::label, Foam::label>>   boundaryCell_;		

	bool 						listConstructed_;

	volScalarField 			boundaryPatchNum_;
	
	/// Members

	//
	bool performCellNeighborSearch()const
	{
		if(cMesh_.moving()) return true;
		if(!listConstructed_)return true;
		return false;
	}


	virtual Foam::scalar boundRatio()const = 0;

	/// construct neighbors of cells based on neighborLength 
	bool cellNeighborsSearch();

	 
public:

	/// Type info
	TypeInfo("statistical");

	/// Construc from dictionary 
	statistical(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		Plus::centerMassField& 	centerMass, 
		Plus::realProcCMField& 	parDiam);

	/// Destructor
	~statistical() override = default ;

	inline
	Foam::scalar neighborLength()const
	{
		return neighborLength_;
	}
	
}; 

} // pFlow::coupling


#endif
