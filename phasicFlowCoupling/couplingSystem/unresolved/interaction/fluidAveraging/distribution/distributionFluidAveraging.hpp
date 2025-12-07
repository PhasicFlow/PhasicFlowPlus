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
 * @class distributionFluidAveraging
 * @brief Averages fluid properties using particle size distribution
 *        weighting.
 * 
 * This concrete implementation of @ref fluidAveraging uses particle
 * distribution functions to compute weighted averages of fluid properties.
 * Instead of using a single cell value, this method distributes each
 * particle's influence across neighboring cells based on Gaussian-like
 * distribution, then performs an inverse distribution to map the result
 * back to particle locations.
 *  
 * The calculate() method uses the coupling system's distribution object
 * to perform the weighted averaging via inverseDistributeValue().
 * 
 * @note Requires cell distribution
 *       (requireCellDistribution() returns true)
 * @see fluidAveraging, cellFluidAveraging, interpolateFluidAveraging
 */

#ifndef __distributionFluidAveraging_hpp__
#define __distributionFluidAveraging_hpp__

#include "fluidAveraging.hpp"

namespace pFlow::coupling
{

class distributionFluidAveraging
:
    public fluidAveraging
{
public:

    TypeInfo("distribution");

    /// Constructor with type, coupling system, and name
    distributionFluidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name);

    /// Destructor
    virtual ~distributionFluidAveraging() = default;
    
    /// Virtual constructor macro for runtime selection
    add_vCtor
    (
        fluidAveraging,
        distributionFluidAveraging,
        word
    );
     
    /// Averages fluid field using particle size distribution weighting
    void calculate(const Foam::volVectorField& orgField);
 
    /// Returns true as distribution averaging requires cell distribution
    bool requireCellDistribution()const override
    {
        return true;
    }
    
};

}


#endif //__distributionFluidAveraging_hpp__
