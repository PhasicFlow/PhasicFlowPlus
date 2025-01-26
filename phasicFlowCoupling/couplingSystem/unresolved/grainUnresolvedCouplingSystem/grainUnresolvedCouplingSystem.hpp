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

#ifndef __grainUnresolvedCouplingSystem_hpp__
#define __grainUnresolvedCouplingSystem_hpp__


#include "UnresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "drag.hpp"

namespace pFlow::coupling
{

template<typename DistributorType>
class grainUnresolvedCouplingSystem
:
	public UnresolvedCouplingSystem<DistributorType>
{
private:

	uniquePtr<porosity>					porosity_ = nullptr;

	uniquePtr<drag<DistributorType>> 	drag_ = nullptr;

	//uniquePtr<velocityInterpolate<interpolate>> interpolate_;

	bool 		requiresDistribution_ = false;

	Plus::realProcCMField		courseGrainFactor_;

	Timer 		porosityTimer_;

	Timer 		interactionTimer_;

protected:

	bool distributeParticleFields() override;

	Plus::realProcCMField& courseGrainFactor()
	{
		return courseGrainFactor_;
	}

public:

	TypeInfoTemplate11("grainUnresolvedCouplingSystem", DistributorType);

	grainUnresolvedCouplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	grainUnresolvedCouplingSystem(const grainUnresolvedCouplingSystem&) = delete;
	
	grainUnresolvedCouplingSystem& operator=(const grainUnresolvedCouplingSystem&) = delete;

	grainUnresolvedCouplingSystem(grainUnresolvedCouplingSystem&&) = delete;

	grainUnresolvedCouplingSystem& operator=(grainUnresolvedCouplingSystem&&) = delete;

	~grainUnresolvedCouplingSystem() override = default;

	add_vCtor
	(
		unresolvedCouplingSystem,
		grainUnresolvedCouplingSystem,
		word
	);

	const Plus::realProcCMField& courseGrainFactor()const
	{
		return courseGrainFactor_;
	}

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
		return "grain";
	}
	
	
};


} // pFlow::coupling


#include "grainUnresolvedCouplingSystem.C"

#endif //__grainUnresolvedCouplingSystem_hpp__
