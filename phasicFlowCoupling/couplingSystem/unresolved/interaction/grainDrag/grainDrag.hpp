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

#ifndef __grainDrag_hpp__ 
#define __grainDrag_hpp__

#include "drag.hpp"
#include "grainUnresolvedCouplingSystem.hpp"
#include "fluidVelocity.hpp"
#include "solidVelocity.hpp"
#include "schedule.hpp"

namespace pFlow::coupling
{

template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
class grainDrag
:
	public drag<DistributorType>
{
public:
	
	using DragType = drag<DistributorType>;

	using GrainDragType = grainDrag<DistributorType, DragClosureType, useCellDistribution>;

private:

	DragClosureType 				dragClosure_;

	const Plus::realProcCMField& 	courseGrainFactor_;

	uniquePtr<fluidVelocity> 	fluidVelocityPtr_ = nullptr;

	uniquePtr<solidVelocity>    solidVelocityPtr_ = nullptr;
	
	Foam::word 			fVelocityType_;
	
	Foam::word 			sVelocityType_;

	void calculateDragForce
	(
		const DragClosureType& 			dragClosure,
		const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce
	);

public:

	// type info
	TypeInfoTemplate211(
		"grainDrag",
		DistributorType::TYPENAME(),
		DragClosureType,
		useCellDistribution?"withDistribution":"noDistribution");

	grainDrag(
		const UnresolvedCouplingSystem<DistributorType>& uCS, 
		const porosity& 				prsty);

	~grainDrag() override = default ;

	add_vCtor
	(
		DragType,
		GrainDragType,
		couplingSystem	
	);

	 
	void calculateDragForce(
		const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce)override;

	bool requireCellDistribution()const override
	{
		if constexpr (useCellDistribution)
			return this->cellDistribution().requireCellDistribution();
		else
			return sVelocityType_ == "cell";
	}
	
	
}; 

} // pFlow::coupling


#include "grainDrag.C"

#endif // __grainDrag_hpp__
