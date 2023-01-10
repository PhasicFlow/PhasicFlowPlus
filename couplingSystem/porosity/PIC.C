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

// from OpenFOAM
#include "fvCFD.H"


#include "PIC.hpp"



pFlow::coupling::PIC::PIC(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam)
{

}


bool pFlow::coupling::PIC::calculatePorosity()
{
	
	auto solidVol = Foam::scalarField(
		this->size(), 
		static_cast<Foam::scalar>(0));


	for(size_t i=0; i<centerMass_.size(); i++)
	{
		
		auto cellId = cMesh_.findCell(centerMass_[i], parCellIndex_[i]);
		if( cellId >= 0 )
		{
			solidVol[cellId] += 
				static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));
		}
		parCellIndex_[i] = cellId;
	}


	this->field() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()) );

	//alpha.ref() = 1 - solidVol/alpha.mesh().V();

	// we may need to update boundary condditions
	// we also need to check if the old time step value is stored or not.

	return true;
}