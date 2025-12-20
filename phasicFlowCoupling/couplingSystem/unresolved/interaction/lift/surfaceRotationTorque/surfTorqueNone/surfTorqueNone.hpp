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
 * @class surfTorqueNone
 * @brief Null surface rotation torque model - no torque applied.
 *
 * Disables torque computation from surface rotation effects.
 */

#ifndef __surfTorqueNone_hpp__
#define __surfTorqueNone_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorque.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfTorqueNone
:
    public surfaceRotationTorque
{


public:

    TypeInfo("none");

    surfTorqueNone(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfTorqueNone() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        surfTorqueNone,
        couplingSystem
    );

    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) override;

};

}

#endif // __surfTorqueNone_hpp__