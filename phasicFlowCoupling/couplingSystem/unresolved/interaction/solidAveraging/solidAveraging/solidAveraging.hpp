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

/// @file solidAveraging.hpp
/// @brief Base class for particle property averaging over fluid cells.
///
/// @class pFlow::coupling::solidAveraging
/// @ingroup interaction
///
/// @brief Abstract interface for solid (particle) property averaging.
///
/// Maps particle-centered properties (velocity, mass, etc)
/// from Lagrangian particles onto the Eulerian fluid grid. Required for
/// calculating cell-averaged solid properties used in momentum coupling
/// (drag, lift, etc.).
///
/// ## Implementation Approaches
/// - **particleSolidAveraging**: Direct particle values without smoothing
/// - **distributionSolidAveraging**: Kernel-based smoothing with back-
///   interpolation to particles
///
/// ## Key Methods
/// - calculate(): Mass-weighted (volume-based) averaging
/// - calculateNumberBased(): Number-weighted averaging
/// - requireCellDistribution(): Whether distribution machinery is used
///
/// @see particleSolidAveraging
/// @see distributionSolidAveraging

#ifndef __solidAveraging_hpp__
#define __solidAveraging_hpp__



// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"
#include "span.hpp"

// from phasicFlowPlus
#include "procCMFields.hpp"
#include "couplingMesh.hpp"
#include "porosity.hpp"


namespace pFlow::coupling
{

class unresolvedCouplingSystem;
class distributionBase;
class porosity;

class solidAveraging
{
protected:

    /// @brief Averaged particle property field reference.
    const Plus::realx3ProcCMField*  parAvFieldRef_;

    /// @brief Reference to the unresolved coupling system.
    const unresolvedCouplingSystem& uCS_;

    /// @brief Reference to the porosity field with particle/cell data.
    const porosity&                 porosity_;

public:

    // type info
	TypeInfo("solidAveraging");

    /// @brief Constructor initializing the averaging system.
    solidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty,
        const word&                     name);

    virtual 
    ~solidAveraging()=default;

    /// @brief Virtual constructor for dynamic object creation.
    create_vCtor
	(
		solidAveraging,
		word,
		(
			const word&                     type,
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty,
            const word&                     name
		),
		(type, uCS, prsty, name)
	);


    /// @brief Perform mass-weighted (volume-based) averaging of particle property.
    virtual 
    void calculate(
        const Plus::realx3ProcCMField& particleField) = 0;
    
    /// @brief Perform number-weighted averaging of particle property.
    virtual 
    void calculateNumberBased(
        const Plus::realx3ProcCMField& particleField) = 0;
    
    /// @brief Check if cell distribution machinery is required.
    virtual 
    bool requireCellDistribution()const = 0;

    /// @brief Get particle property value by particle index.
    inline
    const realx3& value(Foam::scalar celli, size_t parIdx)const
    {
        return (*parAvFieldRef_)[parIdx];
    }

    /// @brief Get the entire averaged property field.
    inline
    const Plus::realx3ProcCMField& field()const
    {
        return *parAvFieldRef_;
    }

    /// @brief Get field data as a span for efficient access.
    inline pFlow::span<realx3> fieldSpan()const
    {
        return pFlow::span<realx3>(const_cast<realx3*>(parAvFieldRef_->data()), parAvFieldRef_->size());
    }

    /// @brief Factory method to create appropriate averaging instance.
    static
    uniquePtr<solidAveraging> create
    (
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty,
        const word&                     name
    );

};

}


#endif //__solidAveraging_hpp__


