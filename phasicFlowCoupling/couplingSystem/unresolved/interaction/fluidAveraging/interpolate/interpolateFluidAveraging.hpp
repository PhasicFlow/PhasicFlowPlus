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
 * @class interpolateFluidAveraging
 * @brief Interpolates fluid properties to particle center of mass
 *        locations.
 * 
 * This concrete implementation of @ref fluidAveraging uses spatial
 * interpolation to compute fluid field values at the exact center.  
 * 
 * The calculate() method performs spatial interpolation using the mesh
 * infrastructure to find neighboring cell values and interpolate at the
 * particle's actual position. 
 * 
 * @note Does not require cell distribution
 *       (requireCellDistribution() returns false)
 * @see fluidAveraging, cellFluidAveraging, distributionFluidAveraging
 */

#ifndef __interpolateFluidAveraging_hpp__
#define __interpolateFluidAveraging_hpp__

#include "fluidAveraging.hpp"

namespace pFlow::coupling
{

class interpolateFluidAveraging
:
    public fluidAveraging
{
public:

    TypeInfo("interpolate");

    /// Constructor with type, coupling system, and name
    interpolateFluidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name);

    /// Destructor
    virtual ~interpolateFluidAveraging() = default;
    
    /// Virtual constructor macro for runtime selection
    add_vCtor
    (
        fluidAveraging,
        interpolateFluidAveraging,
        word
    );
     
    /// Interpolates fluid field to particle center of mass locations
    void calculate(const Foam::volVectorField& orgField);
 
    /// Returns false as interpolation does not require distribution
    bool requireCellDistribution()const override
    {
        return false;
    }
    
};

}


#endif //__interpolateFluidAveraging_hpp__
