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
 * @class Saffmann
 * @brief Low Reynolds number lift model (Re_p < 1).
 *
 * Classic Saffman (1965) model for shear flow lift.
 * @f[ C_{L\omega} = \frac{18}{\pi^2}\sqrt{\frac{Sr}{Re_p}} J -
 *     \frac{11}{8} Sr \exp(-0.5 Re_p) @f]
 * @f[ C_{L\Omega} = Rr @f]
 *
 * @reference Saffman, P.G.T., 1965. J. Fluid Mech. 22, 385–400
 */

#ifndef __Saffmann_hpp__
#define __Saffmann_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "lift.hpp"


namespace pFlow::coupling
{


class Saffmann
:
    public lift
{
private:

    /// Temporary storage for computed lift force field
    tmp<Foam::volVectorField>   tmpLiftForce_;

    /// Minimum Reynolds number threshold to avoid numerical issues
    Foam::scalar        residualRe_;
    
public:

    // type info
    TypeInfo("Saffmann");

    /// Constructs Saffmann lift model from coupling system and porosity
    Saffmann(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~Saffmann() = default;
    
    add_vCtor
    (
        lift,
        Saffmann,
        couplingSystem
    );

    /// Computes lift force using Saffman (1965) correlation
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

#endif // __Saffmann_hpp__