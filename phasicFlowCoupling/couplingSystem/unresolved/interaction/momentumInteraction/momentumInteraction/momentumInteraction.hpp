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

#ifndef __momentumInteraction_hpp__
#define __momentumInteraction_hpp__

/**
 * @class momentumInteraction
 * @brief Manages momentum exchange between fluid and particles in unresolved coupling.
 * 
 * This class handles the calculation and application of momentum coupling between 
 * the fluid phase and particle phase in an unresolved CFD-DEM coupling framework. 
 * It integrates multiple interaction models including drag forces, lift forces, 
 * fluid velocity averaging, and solid phase averaging.
 * 
 * The class supports two momentum exchange strategies:
 * - **Distribution Method**: Distributes momentum exchanges to cells based on 
 *   particle distributions within cells
 * - **Cell Method**: Applies momentum exchanges directly to cell centers
 * 
 * The main responsibilities include:
 * - Computing drag forces on particles via @ref drag model closures (Beetstra, 
 *   DiFelice, ErgunWenYu, etc.)
 * - Computing lift forces on particles via @ref lift models
 * - Averaging fluid properties (velocity) to particle locations
 * - Averaging solid phase properties to cell locations
 * - Managing porosity effects on momentum interactions
 * - Calculating source terms (Su, Sp) for momentum equations
 * 
 * @see drag, lift, fluidAveraging, solidAveraging, porosity
 */

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "Timer.hpp"

// from phasicFlowPlus
#include "procCMFields.hpp"
#include "drag.hpp"
#include "lift.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "PCM.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;
class porosity;


class momentumInteraction
{
private:
    
    /// Reference to the porosity field for cells
    const porosity&            porosity_;

    /// Flag indicating if momentum is exchanged using distribution (true) or cell method (false)
    bool                       momentumExchangeDistribute_;

    /// Flag indicating if cell distribution is required for the coupling
    bool                       requireCellDistribution_ = false;

    /// Pointer to the drag force model
    uniquePtr<drag>            drag_;

    /// Pointer to the lift force model
    uniquePtr<lift>            lift_;

    /// Pointer to the fluid property averaging model (e.g., fluid velocity at particles)
    uniquePtr<fluidAveraging>  fluidAveraging_;

    /// Pointer to the solid phase averaging model (e.g., particle properties to cells)
    uniquePtr<solidAveraging>  solidAveraging_;

    /// Pointer to the particle classification model for non-distributed momentum exchange
    uniquePtr<PCM>             noDistribution_ = nullptr;

    /// Timer for tracking performance of momentum interaction calculations
    Timer                      momentumInteractionTimer_;

   
public:

    TypeInfo("momentumInteraction");

    /// Constructor initializing momentum interaction with unresolved coupling system and porosity field
    momentumInteraction(
        const unresolvedCouplingSystem& uCS,
        const porosity& prsty);

    /// Destructor
    virtual ~momentumInteraction() = default;

    /// Returns the source term (Su) for momentum equation from drag forces
    inline
    const Foam::volVectorField& Su()const
    {
        return std::as_const<const drag&>(*drag_).Su();
    }

    /// Returns the implicit source coefficient (Sp) for momentum equation from drag forces
    inline
    const Foam::volScalarField& Sp()const
    {
        return std::as_const<const drag&>(*drag_).Sp();
    }

    /// Returns the lift force field acting on particles
    inline
    Foam::tmp<Foam::volVectorField> liftForce()const
    {
        return lift_->liftForce();
    }

    /// Returns a const reference to the porosity field
    inline
    const porosity& Porosity()const
    {
        return porosity_;
    }

    /// Returns a const reference to the unresolved coupling system
    const unresolvedCouplingSystem& uCS()const;

    /// Returns true if cell distribution is required for momentum exchange
    inline
    bool requireCellDistribution()const
    {
        return requireCellDistribution_;
    }

    /// Calculates momentum coupling forces and torques from fluid on particles
    virtual
    void calculateCoupling(
        const Foam::volVectorField&     U,
        const Plus::realx3ProcCMField&  vp,
        const Plus::realx3ProcCMField&  wp,
        Plus::realx3ProcCMField&        fluidForce,
        Plus::realx3ProcCMField&        fluidTorque);

    /// Returns a const reference to the momentum interaction dictionary
    const Foam::dictionary& dict()const;

    /// Returns the momentum interaction dictionary from the unresolved coupling system
    static 
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);
};

}

#endif //__momentumInteraction_hpp__