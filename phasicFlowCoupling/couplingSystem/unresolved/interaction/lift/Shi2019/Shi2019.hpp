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
 * @class Shi2019
 * @brief Advanced lift model for wide Reynolds number range.
 *
 * Shi and Rzehak (2019) model with improved accuracy for shear and
 * spin lift across wide Reynolds numbers.
 *
 * **Spin Lift:**
 * @f[ C_{L\Omega} = Rr \left( 1 - 0.62 \tanh(0.3 \sqrt{Re_p}) - 0.24
 *     \frac{\tanh(0.01 Re_p)}{\tanh(0.8 \sqrt{Rr})} \arctan[0.47(Rr-1)]
 *     \right) @f]
 *
 * **Shear Lift:** Re-dependent formulation with improved prediction
 * at high Reynolds numbers.
 *
 * @reference Shi, P., Rzehak, R., 2019. Chem. Eng. Sci. 208, 115145
 */

#ifndef __Shi2019_hpp__
#define __Shi2019_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "lift.hpp"


namespace pFlow::coupling
{


class Shi2019
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
    TypeInfo("Shi2019");

    /// Constructs Shi2019 lift model from coupling system and porosity
    Shi2019(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);
    

    virtual ~Shi2019() = default;
    
    add_vCtor
    (
        lift,
        Shi2019,
        couplingSystem
    );

    /// Computes lift force using Shi and Rzehak (2019) correlation
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

#endif // __Shi2019_hpp__