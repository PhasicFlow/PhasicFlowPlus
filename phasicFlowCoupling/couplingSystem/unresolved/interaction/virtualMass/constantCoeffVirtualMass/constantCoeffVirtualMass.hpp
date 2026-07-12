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
/*
 * @class constantCoeffVirtualMass
 * @brief Virtual (added) mass force model with constant coefficient
 *
 * Models the force required to accelerate the fluid surrounding a particle.
 * This force becomes significant when the fluid density is comparable to 
 * particle density (insignificant when rho_f/rho_p ~ 10^(-3)).
 * The Cvm_ is typically 0.5 for spheres.
 * The virtual mass force is given by:
 * @f[ \mathbf{F}_{vm} = C_{vm} \rho_f V_p \left( \frac{D\mathbf{u}_f}{Dt} - \frac{d\mathbf{v}_p}{dt} \right) @f]

*/

#ifndef __constantCoeffVirtualMass_hpp__
#define __constantCoeffVirtualMass_hpp__

#include "virtualMass.hpp"

namespace pFlow::coupling
{


class constantCoeffVirtualMass 
: 
    public virtualMass
{
private:

    /// Temporary storage for computed virtual mass force field
    tmp<Foam::volVectorField> tmpVirtualMassForce_;
    
    /// constant coefficient of virtual mass
    Foam::scalar 	Cvm_;
    
    /// is flow compressible 
    bool isCompressible_ = false;
        
public:

    TypeInfo("constantCoeff");

    constantCoeffVirtualMass(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~constantCoeffVirtualMass() = default;

    add_vCtor
    (
        virtualMass,
        constantCoeffVirtualMass,
        couplingSystem
    );

    void calculateVirtualMassForce
    (
        const Foam::volVectorField& U,
        const Plus::realx3ProcCMField& parAcc,
        const Plus::realProcCMField& diameter,
        Plus::realx3ProcCMField& particleForce
    ) override;
    
    /// Returns the computed virtual mass force field
    Foam::tmp<Foam::volVectorField> virtualMassForce()const override
    {
        return tmpVirtualMassForce_;
    }
    
    inline 
    bool isCompressible()const
    {
        return isCompressible_;
    }
    
};

} // pFlow::coupling

#endif
