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

#ifndef __powerDistributionKernel_hpp__ 
#define __powerDistributionKernel_hpp__


// from phasicFlow-coupling
#include "statistical.hpp"


namespace pFlow::coupling
{

/**
 * powerDistributionKernel model for porosity calculation
 * 
 * This model distributes particle volume to the target cell and neighboring 
 * cell based on the a Gaussian distribution and the distance between 
 * particle center and cell center. 
 * 
 */
class powerDistributionKernel
: 
	public statistical
{
protected:

	Foam::Switch 				filterEmpty_;

	Foam::scalar boundRatio()const override
	{
		return 1.0;
	} 		

public:

	/// Type info
	TypeInfo("powerDistributionKernel");

	/// Construc from dictionary 
	powerDistributionKernel(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		Plus::centerMassField& 	centerMass, 
		Plus::realProcCMField& 	parDiam);

	/// Destructor
	virtual ~powerDistributionKernel() = default;

	/// Add this constructor to the list of virtual constructors
	add_vCtor
	(
		porosity,
		powerDistributionKernel,
		dictionary
	);

	bool internalFieldUpdate() override;


}; 

} // pFlow::coupling


#endif
