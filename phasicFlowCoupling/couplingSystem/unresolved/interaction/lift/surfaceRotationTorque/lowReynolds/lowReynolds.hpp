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
 * @class lowReynolds
 * @brief Low Reynolds number surface rotation torque model.
 *
 * Happel & Brenner (1973) model with constant correction factors
 * (@f$ f_{spin} = 1.0 @f$, @f$ f_{shear} = 1.0 @f$) for
 * viscous-dominated flows.
 *
 * @reference Happel, J., Brenner, H., 1973. Low Reynolds Number
 * Hydrodynamics. Noordhoff International.
 */

#ifndef __lowReynolds_hpp__
#define __lowReynolds_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

#include "surfaceRotationTorque.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class lowReynolds
:
    public surfaceRotationTorque
{

public:

    TypeInfo("lowReynolds");

    /// Constructs lowReynolds torque model from coupling system and porosity
    lowReynolds(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~lowReynolds() = default;

    add_vCtor
    (
        surfaceRotationTorque,
        lowReynolds,
        couplingSystem
    );

    /// Computes surface torque using Happel & Brenner (1973) model
    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) override;

};

}

#endif // __lowReynolds_hpp__