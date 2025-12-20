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
 * @class none
 * @brief Null lift force model - no lift forces applied.
 *
 * Disables lift force computation. Use when lift effects are
 * negligible.
 */

#ifndef __none_hpp__
#define __none_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "lift.hpp"


namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class none
:
    public lift
{
private:

    /// Temporary storage for lift force field (always zero)
    tmp<Foam::volVectorField>   tmpLiftForce_;
    
public:

    // type info
    TypeInfo("none");

    /// Constructs null lift model from coupling system and porosity
    none(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~none() = default;
    
    add_vCtor
    (
        lift,
        none,
        couplingSystem
    );

    /// No-op method - does not compute lift forces
    void calculateLiftForce(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque) override;

    /// Returns zero lift force field
    Foam::tmp<Foam::volVectorField> liftForce()const override
    {
        return tmpLiftForce_;
    }

};
    
}

#endif // __none_hpp__