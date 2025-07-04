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
#ifndef __GaussianIntegral_hpp__ 
#define __GaussianIntegral_hpp__


#include "distribution.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class GaussianIntegral
:
	public distribution
{
private:
	
	// radius of circule for cell neighbor list 
	Foam::label 				maxLayers_;

public:

	/// Type info
	TypeInfoNV("GaussianIntegral");

	/// Construc from dictionary 
	GaussianIntegral(
		Foam::dictionary 		dict, 
		const couplingMesh& 	cMesh,
		const Plus::centerMassField& centerMass);

	/// Destructor
	~GaussianIntegral() = default;

	void updateWeights
	(
		const Plus::procCMField<Foam::label> & parCellIndex,
		const Plus::procCMField<real> & parDiameter
	);

	bool requireCellDistribution()const 
	{
		return true;
	}

}; 

} // pFlow::coupling


#endif
