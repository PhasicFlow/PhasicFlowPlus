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

/**
 * @class sphereDrag
 * @brief Templated drag model for sphere-scale unresolved coupling.
 *
 * This template class implements drag force calculations for spherical
 * particles in unresolved Eulerian-Lagrangian coupling, where the drag
 * closure is selected at compile time via template parameter.
 *
 * @tparam DragClosureType The drag closure model (e.g., DiFelice,
 *                         ErgunWenYu, Beetstra, Rong).
 *
 * @details
 * Computes drag forces for spheres by applying the selected drag closure
 * model to calculate particle-fluid momentum exchange. Inherits from the
 * drag base class and is used in sphere-scale unresolved simulations.
 */

#ifndef __sphereDrag_hpp__ 
#define __sphereDrag_hpp__

#include "drag.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "distributionBase.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

template<typename DragClosureType>
class sphereDrag
:
	public drag
{
public:

    /// @brief Template type alias for sphere drag implementation.
	using SphereDragType = sphereDrag<DragClosureType>;

private:

    /// @brief Drag closure model instance.
    DragClosureType 			dragClosure_;

public:

    // type info
    TypeInfoTemplate11("sphereDrag",DragClosureType);

    /// @brief Constructor from unresolved coupling system.
    sphereDrag(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);

    /// @brief Destructor.
    virtual ~sphereDrag() override = default ;

    add_vCtor
    (
        drag,
        SphereDragType,
        couplingSystem	
    );

    /// @brief Calculate drag force for spherical particles.
    void calculateDragForce(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,    
        Plus::realx3ProcCMField& 		particleForce)override;

}; 

} // pFlow::coupling


#include "sphereDrag.C"

#endif // __sphereDrag_hpp__
