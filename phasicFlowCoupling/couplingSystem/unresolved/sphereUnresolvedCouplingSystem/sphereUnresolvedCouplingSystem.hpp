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

#ifndef __sphereUnresolvedCouplingSystem_hpp__
#define __sphereUnresolvedCouplingSystem_hpp__


#include "UnresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "drag.hpp"

namespace pFlow::coupling
{

template<typename DistributorType>
class sphereUnresolvedCouplingSystem
:
	public UnresolvedCouplingSystem<DistributorType>
{
private:

	uniquePtr<porosity>								porosity_ = nullptr;

	uniquePtr<drag<DistributorType>> 	drag_ = nullptr;

	//uniquePtr<velocityInterpolate<interpolate>> interpolate_;

	bool 		requiresDistribution_ = false;

	Timer 		porosityTimer_;

	Timer 		interactionTimer_;

public:

	TypeInfoTemplate11("sphereUnresolvedCouplingSystem", DistributorType);

	sphereUnresolvedCouplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	sphereUnresolvedCouplingSystem(const sphereUnresolvedCouplingSystem&) = delete;
	
	sphereUnresolvedCouplingSystem& operator=(const sphereUnresolvedCouplingSystem&) = delete;

	sphereUnresolvedCouplingSystem(sphereUnresolvedCouplingSystem&&) = delete;

	sphereUnresolvedCouplingSystem& operator=(sphereUnresolvedCouplingSystem&&) = delete;

	~sphereUnresolvedCouplingSystem() override = default;

	add_vCtor
	(
		unresolvedCouplingSystem,
		sphereUnresolvedCouplingSystem,
		word
	);


	void calculateFluidInteraction() override;

	void calculatePorosity() override;

	const Foam::volScalarField& alpha()const override
	{
		return porosity_().alpha();
	}

	const Foam::volScalarField& Sp()const override
	{
		return drag_().Sp();
	}
	
	const Foam::volVectorField& Su()const override
	{
		return drag_().Su();
	}

	word shapeTypeName() const override
	{
		return "sphere";
	}
	
	
};


} // pFlow::coupling


#include "sphereUnresolvedCouplingSystem.C"

#endif //__sphereUnresolvedCouplingSystem_hpp__

