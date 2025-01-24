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

#ifndef __sphereDrag_hpp__ 
#define __sphereDrag_hpp__

#include "drag.hpp"


namespace pFlow::coupling
{

template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
class sphereDrag
:
	public drag<DistributorType>
{
public:
	
	using DragType = drag<DistributorType>;

	using SphereDragType = sphereDrag<DistributorType, DragClosureType, useCellDistribution>;

private:

	DragClosureType 			dragClosure_;


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
		"sphereDrag",
		DistributorType::TYPENAME(),
		DragClosureType,
		useCellDistribution?"withDistribution":"noDistribution");

	sphereDrag(
		const UnresolvedCouplingSystem<DistributorType>& uCS, 
		const porosity& 				prsty);

	~sphereDrag() override = default ;

	add_vCtor
	(
		DragType,
		SphereDragType,
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
			return true;
		else
			return false;
	}
	
	
}; 

} // pFlow::coupling


#include "sphereDrag.C"

#endif // __sphereDrag_hpp__
