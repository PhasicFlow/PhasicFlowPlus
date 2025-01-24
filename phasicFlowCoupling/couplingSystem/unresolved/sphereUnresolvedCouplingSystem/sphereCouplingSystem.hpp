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

#ifndef __sphereCouplingSystem_hpp__
#define __sphereCouplingSystem_hpp__


#include "couplingSystem.hpp"
#include "porosity.hpp"
#include "drag.hpp"

namespace pFlow::coupling
{


class sphereCouplingSystem
:
	public couplingSystem
{
private:
	
	uniquePtr<porosity>			porosity_ = nullptr;

	uniquePtr<drag> 			drag_ = nullptr;

	Timer 						porosityTimer_;

	Timer 						interactionTimer_;


public:


	sphereCouplingSystem(
		word demSystemName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	sphereCouplingSystem(const sphereCouplingSystem&) = delete;
	
	sphereCouplingSystem& operator=(const sphereCouplingSystem&) = delete;

	sphereCouplingSystem(sphereCouplingSystem&&) = delete;

	sphereCouplingSystem& operator=(sphereCouplingSystem&&) = delete;

	~sphereCouplingSystem() override = default;

	void calculateFluidInteraction() override;

	void calculatePorosity() override;

	const Foam::volScalarField& alpha()const override;

	const Foam::volScalarField& Sp()const override;
	
	const Foam::volVectorField& Su()const override;


}; 

} // pFlow::coupling

#endif //__sphereCouplingSystem_hpp__
