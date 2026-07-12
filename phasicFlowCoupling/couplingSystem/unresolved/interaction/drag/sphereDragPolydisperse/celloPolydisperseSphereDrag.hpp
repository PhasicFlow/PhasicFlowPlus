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
 * @class celloPolydisperseSphereDrag
 * @brief Drag model for cell-based polydisperse sphere simulations.
 *
 */

#ifndef __celloPolydisperseSphereDrag_hpp__ 
#define __celloPolydisperseSphereDrag_hpp__

#include "polydisperseDragBase.hpp"

namespace pFlow::coupling
{

class celloPolydisperseSphereDrag
:
    public polydisperseDragBase
{

public:

    // type info
    TypeInfo("sphereDrag<celloPolydisperse>");

    /// @brief Constructor from unresolved coupling system.
    celloPolydisperseSphereDrag(
        const unresolvedCouplingSystem& uCS, 
        const porosity&                 prsty);

    /// @brief Destructor.
    virtual ~celloPolydisperseSphereDrag() override = default;
    
    add_vCtor
    (
        drag,
        celloPolydisperseSphereDrag,
        couplingSystem
    )

    /// @brief Calculate drag force for polydisperse spheres.
    virtual 
    void calculateDragForce(
        const fluidAveraging&           fluidVelocity,
        const solidAveraging&           parVelocity,
        const Plus::realProcCMField&    diameter,
        const distributionBase&         cellDistribution,    
        Plus::realx3ProcCMField&        particleForce,
        Foam::volScalarField&           Sp,
        Foam::volVectorField&           Su) override;
};

}

#endif // __celloPolydisperseSphereDrag_hpp__