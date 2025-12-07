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
 * @class subDivision9Dist
 * @brief 9-sub-volume distribution method for particle-to-cell coupling.
 *
 * Divides each particle sphere into 9 equal volumetric segments and locates 
 * each segment in the appropriate cell. Provides good balance between accuracy 
 * and computational efficiency for porosity calculations in unresolved simulations.
 */

#ifndef __subDivision9Dist_hpp__
#define __subDivision9Dist_hpp__

#include "distribution.hpp"

namespace pFlow::coupling
{

class couplingMesh;

class subDivision9Dist
:
    public distribution
{

public:

    /// Type info
    TypeInfo("subDivision9");

    /// Initialize 9-sub-volume distribution from dictionary
    subDivision9Dist(
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~subDivision9Dist() = default;

    add_vCtor
    (
        distributionBase,
        subDivision9Dist,
        dictionary
    );

    /// Update distribution weights by subdividing particles into 9 volumes
    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "subDivision9";
    }

    

}; // end subDivision9Dist

} // end pFlow::coupling


#endif //__subDivision9Dist_hpp__
