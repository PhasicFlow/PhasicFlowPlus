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
 * @class Loth2008
 * @brief Wide-range lift model (Re_p < 5000).
 *
 * Loth (2008) model accounting for shear and spin effects with
 * negative lift prediction capability.
 *
 * **Spin Lift:**
 * @f[ C_{L\Omega} = Rr \left( 1 - \{0.675 +
 *     0.15(1 + \tanh[0.28(Rr-2)])\} \tanh(0.18 \sqrt{Re_p})
 *     \right) @f]
 *
 * **Shear Lift (Re_p ≤ 50):** Complex function via @f$ J(\epsilon) @f$
 * **Shear Lift (Re_p > 50):** Negative lift correlation
 *
 * @reference Loth, E., 2008. AIAA Journal, 46(4), 801–809
 */

#ifndef __Loth2008_hpp__
#define __Loth2008_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "lift.hpp"


namespace pFlow::coupling
{


class Loth2008
:
    public lift
{
private:

    /// Temporary storage for computed lift force field
    tmp<Foam::volVectorField>   tmpLiftForce_;

    /// Minimum Reynolds number threshold to avoid numerical issues
    Foam::scalar                residualRe_;
    
public:

    // type info
    TypeInfo("Loth2008");

    /// Constructs Loth2008 lift model from coupling system and porosity
    Loth2008(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~Loth2008() = default;
    
    add_vCtor
    (
        lift,
        Loth2008,
        couplingSystem
    );

    /// Computes lift force using Loth2008 correlation
    void calculateLiftForce(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque) override;

    /// Returns the computed lift force field
    Foam::tmp<Foam::volVectorField> liftForce()const override
    {
        return tmpLiftForce_;
    }

};
    
}

#endif // __Loth2008_hpp__