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
 * @class fluidAveraging
 * @brief Abstract base class for averaging fluid field properties to
 *        particle locations.
 * 
 * This class provides an interface for computing and storing interpolated
 * or averaged fluid properties (such as velocity) at particle center. 
 * Different derived classes implement various strategies
 * for averaging fluid fields:
 * - Cell-based: Uses value from the cell containing the particle
 * - Distribution-based: Distributes cell values weighted by particle
 *   distributions
 * - Interpolation-based: Interpolates field values to exact particle
 *   positions
 * 
 * The averaged field is stored in @ref averagedField_ as a per-particle
 * processor field, which can be accessed via value(), field(), fieldSpan() 
 * methods.
 * 
 * @see cellFluidAveraging, distributionFluidAveraging,
 *      interpolateFluidAveraging
 */

#ifndef __fluidAveraging_hpp__
#define __fluidAveraging_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

#include "procCMFields.hpp"
#include "virtualConstructor.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;


class fluidAveraging
{
protected:
    
    /// Interpolated fluid field at particle's center of mass
    pFlow::Plus::procCMField<Foam::vector>   averagedField_;
    
    /// Reference to unresolved coupling system
    const unresolvedCouplingSystem&          uCS_;

public:

    /// Type information for runtime selection
    TypeInfo("fluidAveraging");

    /// Constructor with type, coupling system, and name
    fluidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name);

    /// Destructor
    virtual ~fluidAveraging() = default;

    /// Virtual constructor macro for runtime selection
    create_vCtor
    (
        fluidAveraging,
        word,
        (
            const word&                     type,
            const unresolvedCouplingSystem& uCS,
            const word&                     name
        ),
        (type, uCS, name)
    );
    
    /// Calculates averaged field from the original field
    virtual 
    void calculate(const Foam::volVectorField& orgField)=0;

    /// Returns the averaged fluid property at particle location
    inline
    const Foam::vector& value(Foam::scalar celli, size_t parIdx)const
    {
        return averagedField_[parIdx];;
    }

    /// Returns the entire averaged field
    inline
    const pFlow::Plus::procCMField<Foam::vector>& field()const
    {
        return averagedField_;
    }

    /// Returns the averaged field as a span for efficient access
    inline pFlow::span<Foam::vector> fieldSpan()const
    {
        return pFlow::span<Foam::vector>(const_cast<Foam::vector*>(averagedField_.data()), averagedField_.size());
    }

    /// Returns true if cell distribution is required for averaging
    virtual 
    bool requireCellDistribution()const = 0;

    /// Creates and returns a new averaging model instance
    static 
    uniquePtr<fluidAveraging> create(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name
    );
    
};

}


#endif //__fluidAveraging_hpp__
