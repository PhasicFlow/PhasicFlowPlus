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
 * @class lift
 * @brief Base class for lift force models.
 *
 * Abstract interface for computing particle lift forces from
 * fluid vorticity and rotation effects.
 * @see Saffmann, Loth2008, Shi2019
 */

#ifndef __lift_hpp__
#define __lift_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from PhasicFlowPlus
#include "procCMFields.hpp"
#include "porosity.hpp"
#include "surfaceRotationTorque.hpp"


namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class lift
{
private:
    
    /// Reference to porosity object for momentum coupling
    const porosity&             porosity_;

    /// Flag to enable/disable lift force field to the output file
    Foam::Switch                printLift_;
    
    /// Flag indicating whether lift force calculation is active
    bool                        liftActive_ = false;

protected:

    /// Surface rotation torque model for particle angular acceleration
    uniquePtr<surfaceRotationTorque> surfaceRotationTorque_;
    

    /// Sets the active status of lift force calculation
    void setLiftActive(bool active)
    {
        liftActive_ = active;
    }

    /// Returns true if lift force output is enabled
    bool printLift()const
    {
        return printLift_;
    }
public:

    // type info
    TypeInfo("lift");

    /// Constructs lift model from unresolved coupling system and porosity
    lift(
        const unresolvedCouplingSystem& uCS, 
        const porosity&                 prsty);
    

    virtual ~lift() = default;

    create_vCtor
    (
        lift,
        couplingSystem,
        (
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty
        ),
        (uCS, prsty)
    );

    /// Calculates lift force and torque, then distributes to fluid cells
    void calculateLiftForceTorque(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque);

    /// Pure virtual method for computing lift force on particles
    virtual 
    void calculateLiftForce(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  parVel,
        const Plus::realx3ProcCMField&  parRotVel,
        const Plus::realProcCMField&    diameter,
        Plus::realx3ProcCMField&        particleForce,
        Plus::realx3ProcCMField&        particleTorque) = 0;
    
    /// Pure virtual method returning lift force field
    virtual 
    Foam::tmp<Foam::volVectorField> liftForce()const =0;
    
    /// Returns true if lift force calculation is active
    inline bool liftActive()const
    {
        return liftActive_;
    }

    /// Returns constant reference to porosity field
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

    /// Returns the configuration dictionary for this model
    const Foam::dictionary& dict()const;
    
    /// Computes fluid vorticity field from velocity
    static
    Foam::tmp<Foam::volVectorField> fluidVorticity(const Foam::volVectorField& U);

    /// Retrieves lift model configuration dictionary
    static
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);

    /// Factory method for creating lift model instances
    static
    uniquePtr<lift> create
    (
        const unresolvedCouplingSystem& uCS, 
        const porosity&                 prsty
    );
};
    
}

#endif // __lift_hpp__