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

/// @file particleSolidAveraging.hpp
/// @brief Direct particle property mapping without spatial smoothing.
///
/// @class pFlow::coupling::particleSolidAveraging
/// @ingroup interaction
/// @extends pFlow::coupling::solidAveraging
///
/// @brief Direct particle property implementation (Particle Centroid Method).
///
/// Returns particle properties directly without spatial averaging or
/// distribution to neighboring cells. Simplest approach with minimal
/// computational cost.
///
/// ## Characteristics
/// - No spatial smoothing
/// - Direct value mapping (identity operation)
/// - Minimal overhead
/// - Does not require distribution machinery
///
/// @see solidAveraging (base class)
/// @see distributionSolidAveraging (with kernel smoothing)

#ifndef __particleSolidAveraging_hpp__
#define __particleSolidAveraging_hpp__

#include "solidAveraging.hpp"

namespace pFlow::coupling
{

class particleSolidAveraging
:
    public solidAveraging
{
public:

    // type info
    TypeInfo("particle");

    /// @brief Constructor for direct particle property averaging.
    particleSolidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty,
        const word&                     name
    );

    /// @brief Destructor.
    ~particleSolidAveraging()override =default;

    /// @brief Virtual constructor macro for factory registration.
    add_vCtor
	(
		solidAveraging,
		particleSolidAveraging,
		word	
	);
    
    /// @brief Perform mass-weighted averaging (returns particle field directly).
    void calculate(const Plus::realx3ProcCMField& particleField) override;
    
    /// @brief Perform number-weighted averaging (returns particle field directly).
    void calculateNumberBased(const Plus::realx3ProcCMField& particleField) override;
    
    /// @brief No distribution machinery required (direct particle method).
    bool requireCellDistribution()const override
    {
        return false;
    }

};

}

#endif //__particleSolidAveraging_hpp__