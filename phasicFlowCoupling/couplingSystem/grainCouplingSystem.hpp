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

#ifndef __grainCouplingSystem_hpp__
#define __grainCouplingSystem_hpp__


#include "couplingSystem.hpp"
#include "porosity.hpp"
#include "grainDrag.hpp"

namespace pFlow::coupling
{


class grainCouplingSystem
:
	public couplingSystem
{
private:
	
	uniquePtr<porosity>			porosity_ = nullptr;

	uniquePtr<grainDrag>		drag_ = nullptr;

	Plus::realProcCMField		courseGrainFactor_;

	Timer 						porosityTimer_;

	Timer 						interactionTimer_;

protected:

	bool distributeParticleFields() override;

	auto& courseGrainFactor()
	{
		return courseGrainFactor_;
	}

public:


	grainCouplingSystem(
		word demSystemName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	grainCouplingSystem(const grainCouplingSystem&) = delete;
	
	grainCouplingSystem& operator=(const grainCouplingSystem&) = delete;

	grainCouplingSystem(grainCouplingSystem&&) = delete;

	grainCouplingSystem& operator=(grainCouplingSystem&&) = delete;

	~grainCouplingSystem() override = default;

	void calculateFluidInteraction() override;

	void calculatePorosity() override;

	const Foam::volScalarField& alpha()const override;

	const Foam::volScalarField& Sp()const override;
	
	const Foam::volVectorField& Su()const override;


}; 

} // pFlow::coupling

#endif //__grainCouplingSystem_hpp__
