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
 * @class cellFluidAveraging
 * @brief Averages fluid properties using cell-center values.
 * 
 * This concrete implementation of @ref fluidAveraging assigns the fluid
 * field value of the cell containing each particle to that particle.
 *  
 * @note Does not require cell distribution
 *       (requireCellDistribution() returns false)
 * @see fluidAveraging, distributionFluidAveraging,
 *      interpolateFluidAveraging
 */

#ifndef __cellFluidAveraging_hpp__
#define __cellFluidAveraging_hpp__

#include "fluidAveraging.hpp"

namespace pFlow::coupling
{

//class couplingMesh;

class cellFluidAveraging
:
    public fluidAveraging
{
public:

    TypeInfo("cell");

    /// Constructor with type, coupling system, and name
    cellFluidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const word&                     name);
    
    /// Destructor
    virtual ~cellFluidAveraging() = default;
    
    /// Virtual constructor macro for runtime selection
    add_vCtor
    (
        fluidAveraging,
        cellFluidAveraging,
        word
    );

    /// Copies fluid field values from cell centers to particles
    void calculate(const Foam::volVectorField& orgField)override;

    /// Returns false as cell averaging does not require distribution
    bool requireCellDistribution()const override
    {
        return false;
    }

};

}


#endif //__cellFluidAveraging_hpp__
