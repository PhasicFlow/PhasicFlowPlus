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
 * @class surfTorqueLoth2008
 * @brief Surface rotation torque model based on Loth (2008).
 *
 * Reynolds-number-dependent correction factors for intermediate to high
 * Reynolds number flows:
 * @f[ f_{spin} = 1.0 + \frac{5}{64\pi} Re_\Omega^{0.6} @f]
 * @f[ f_{shear} = f_{spin} (1 - 0.0075 Re_\omega)
 *     (1 - 0.062\sqrt{Re_p} - 0.001 Re_p) @f]
 *
 * @reference Loth, E., 2008. AIAA Journal, 46(4), 801–809
 */

#ifndef __surfTorqueLoth2008_hpp__
#define __surfTorqueLoth2008_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorque.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfTorqueLoth2008
:
    public surfaceRotationTorque
{

    /// Minimum Reynolds number threshold to avoid numerical issues
    Foam::scalar    residualRe_;

public:

    TypeInfo("Loth2008");

    /// Constructs Loth2008 torque model from coupling system and porosity
    surfTorqueLoth2008(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfTorqueLoth2008() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        surfTorqueLoth2008,
        couplingSystem
    );

    /// Computes surface torque using Loth (2008) correlation
    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) override;

};

}

#endif // __surfTorqueLoth2008_hpp__