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
 * @class subDivision29Dist
 * @brief 29-sub-volume distribution method for particle-to-cell coupling.
 *
 * Divides each particle sphere into 29 equal volumetric segments and locates 
 * each segment in the appropriate cell. Provides high accuracy for porosity 
 * calculations in dense suspensions at the cost of increased computational expense.
 */

#ifndef __subDivision29Dist_hpp__
#define __subDivision29Dist_hpp__

#include "distribution.hpp"

namespace pFlow::coupling
{

class couplingMesh;

class subDivision29Dist
:
    public distribution
{

public:

    /// Type info
    TypeInfo("subDivision29");

    /// Initialize 29-sub-volume distribution from dictionary
    subDivision29Dist(
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~subDivision29Dist() = default;

    add_vCtor
    (
        distributionBase,
        subDivision29Dist,
        dictionary
    );

    /// Update distribution weights by subdividing particles into 29 volumes
    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "subDivision29";
    }

    

}; // end subDivision29Dist

} // end pFlow::coupling


#endif //__subDivision29Dist_hpp__
