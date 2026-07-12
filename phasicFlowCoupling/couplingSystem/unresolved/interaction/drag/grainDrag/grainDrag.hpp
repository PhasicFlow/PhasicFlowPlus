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
 * @class grainDrag
 * @brief Templated drag model for grain-scale unresolved coupling.
 *
 * This template class implements drag force calculations for grain-scale
 * particles in unresolved Eulerian-Lagrangian coupling, where the drag
 * closure is selected at compile time via template parameter.
 *
 * @tparam DragClosureType The drag closure model (e.g., DiFelice,
 *                         ErgunWenYu, Beetstra, Rong).
 *
 */

#ifndef __grainDrag_hpp__ 
#define __grainDrag_hpp__

#include "drag.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "distributionBase.hpp"
#include "momentumGrainUnresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

template<typename DragClosureType>
class grainDrag
:
    public drag
{
public:

    /// @brief Template type alias for grain drag implementation.
    using grainDragType = grainDrag<DragClosureType>;

private:

    /// @brief Drag closure model instance.
    DragClosureType 			dragClosure_;

    /// @brief Coarse-grain factor 
    const Plus::realProcCMField& 	courseGrainFactor_;

public:

    // type info
    TypeInfoTemplate11("grainDrag",DragClosureType);

    /// @brief Constructor from unresolved coupling system.
    grainDrag(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);

    /// @brief Destructor.
    virtual ~grainDrag() override = default ;

    add_vCtor
    (
        drag,
        grainDragType,
        couplingSystem	
    );

    /// @brief Calculate drag force for grain-scale particles.
    void calculateDragForce(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,    
        Plus::realx3ProcCMField& 		particleForce)override;

}; 

} // pFlow::coupling


#include "grainDrag.C"

#endif // __grainDrag_hpp__
