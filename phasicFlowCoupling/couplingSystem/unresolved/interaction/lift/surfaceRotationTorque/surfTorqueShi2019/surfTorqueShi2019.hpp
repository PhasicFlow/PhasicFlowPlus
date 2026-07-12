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
 * @class surfTorqueShi2019
 * @brief Advanced surface rotation torque model based on Shi2019.
 *
 * Improved empirical correlations for wide Reynolds number range:
 * @f[ f_{spin} = 1.0 + \frac{5}{64\pi} Re_\Omega^{0.6} @f]
 * @f[ f_{shear} = f_{spin} [(1 + 0.4 (\exp(-0.0135 Re_\omega) - 1))
 *     (1 - 0.0702 Re_p^{0.455})] @f]
 *
 * Latest model with enhanced accuracy across wide Reynolds number ranges.
 *
 * @reference Shi, P., Rzehak, R., 2019. Chem. Eng. Sci. 208, 115145
 */

#ifndef __surfTorqueShi2019_hpp__
#define __surfTorqueShi2019_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorque.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfTorqueShi2019
:
    public surfaceRotationTorque
{

    /// Minimum Reynolds number threshold to avoid numerical issues
    Foam::scalar    residualRe_;

public:

    TypeInfo("Shi2019");

    /// Constructs Shi2019 torque model from coupling system and porosity
    surfTorqueShi2019(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfTorqueShi2019() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        surfTorqueShi2019,
        couplingSystem
    );

    /// Computes surface torque using Shi and Rzehak (2019) correlation
    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) override;

};

}

#endif // __surfTorqueShi2019_hpp__