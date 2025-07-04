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
#ifndef __Gaussian_hpp__ 
#define __Gaussian_hpp__

// from PhasicFlowPlus
#include "distribution.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class Gaussian
:
	public distribution
{	
private:
	
	// radius of circle for cell neighbor list 
	Foam::scalar 				standardDeviation_;

	Foam::label 				maxLayers_;

public:

	/// Type info
	TypeInfoNV("Gaussian");

	/// Construct from dictionary 
	Gaussian(
		Foam::dictionary 		dict, 
		const couplingMesh& 	cMesh,
		const Plus::centerMassField& centerMass);

	/// Destructor
	~Gaussian() = default;

	inline
	auto standardDeviation()const
	{
		return standardDeviation_;
	}

	void updateWeights(
		const Plus::procCMField<Foam::label> & parCellIndex,
		const Plus::procCMField<real> & parDiameter);
	
	bool requireCellDistribution()const 
	{
		return true;
	}
}; 

} // pFlow::coupling


#endif
