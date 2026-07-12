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
 * @class surfaceRotationTorque
 * @brief Base class for surface rotation torque models.
 *
 * Abstract interface for torque from fluid vorticity and particle
 * rotation interaction. General formula:
 * @f[ \mathbf{T} = -\pi \rho \nu d_p^3 (f_{spin} \boldsymbol{\omega}_p -
 *     0.5 f_{shear} \boldsymbol{\omega}_f) @f]
 */

#ifndef __surfaceRotationTorque_hpp__
#define __surfaceRotationTorque_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "porosity.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class surfaceRotationTorque
{

    /// Reference to porosity object
    const porosity&      porosity_;

public:

    TypeInfo("surfaceRotationTorque");

    /// Constructs surface torque model from coupling system and porosity
    surfaceRotationTorque(
        const unresolvedCouplingSystem& uCS, 
	    const porosity&                 prsty
    );

    virtual ~surfaceRotationTorque() = default;

    create_vCtor
    (
        surfaceRotationTorque,
        couplingSystem,
        (
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty
        ),
        (uCS, prsty)
    );

    /// Returns constant reference to porosity object
    inline
    const porosity& Porosity()const
    {
        return porosity_;
    }

    /// Returns particle-to-cell index mapping
    inline
    const auto& parCellIndex()const
    {
        return porosity_.parCellIndex();
    }

    /// Returns constant reference to the computational mesh
    inline
    const Foam::fvMesh& mesh()const
    {
        return porosity_.mesh();
    }

    /// Returns coupling mesh for particle-cell operations
    inline 
    const auto& cMesh()const
    {
        return porosity_.cMesh();
    }

    /// Returns porosity (fluid volume fraction) field
    const Foam::volScalarField& alpha()const
    {
        return static_cast<const Foam::volScalarField&>(porosity_);
    }

    /// Pure virtual method for computing surface rotation torque
    virtual 
    void calculateSurfaceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleTorque) = 0;

    /// Returns the configuration dictionary for this model
    const Foam::dictionary& dict()const;
    
    /// Retrieves torque model configuration dictionary
    static
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);

    /// Factory method for creating surface torque model instances
    static
    uniquePtr<surfaceRotationTorque> create
    (
        const unresolvedCouplingSystem& uCS, 
        const porosity&                 prsty
    );
};

}

#endif // __surfaceRotationTorque_hpp__