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
 * @class polydisperseDragBase
 * @brief Base class for polydisperse particle drag models.
 *
 * This class provides a foundation for drag force models that handle
 * particles with a distribution of sizes (polydisperse systems).
 *
 */

#ifndef __polydisperseDragBase_hpp__ 
#define __polydisperseDragBase_hpp__

#include "drag.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class polydisperseDragBase
:
	public drag
{
private:

    /// @brief Diameter class index for each particle.
    Plus::int32ProcCMField      diameterClass_;

    /// @brief Flag to control average diameter field output.
    Foam::Switch                writeAverageDiameter_;
    
    /// @brief Volume-averaged diameter field per cell.
    Foam::volScalarField        averageDiameter_;

public:

    // type info
    TypeInfo("polydisperseDragBase");

    /// @brief Constructor from unresolved coupling system.
    polydisperseDragBase(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);

    /// @brief Destructor.
    virtual ~polydisperseDragBase() override = default ;

    /// @brief Get mutable reference to average diameter field.
    Foam::volScalarField& averageDiameter()
    {
        return averageDiameter_;
    }

    /// @brief Get constant reference to average diameter field.
    const Foam::volScalarField& averageDiameter()const
    {
        return averageDiameter_;
    }


}; 

} // pFlow::coupling


#endif // __polydisperseDragBase_hpp__
