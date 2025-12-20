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
 * @class PCM
 * @brief Particle Centroid Method (PCM) for particle-to-cell coupling.
 *
 * The simplest distribution approach with no smoothing. All particle data 
 * is assigned directly to the cell containing the particle center.
 */

#ifndef __PCM_hpp__ 
#define __PCM_hpp__

#include "distributionBase.hpp"

namespace pFlow::coupling
{

class PCM
:
    public distributionBase
{
protected:
    
    /// Construct neighbor lists with specified search length and max layers
    void constructLists(
        const Foam::scalar searchLen, 
        const Foam::label maxLayers) override;
    
    /// Construct neighbor lists with max layers only
    void constructLists(const Foam::label maxLayers)override;

public:    
    
    /// Type info
    TypeInfo("PCM");

    /// Initialize PCM distribution from dictionary
    PCM(
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Initialize PCM distribution without dictionary
    PCM(
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~PCM() = default;

    add_vCtor
    (
        distributionBase,
        PCM,
        dictionary  
    );

    /// Update distribution weights (no smoothing for PCM)
    void updateWeights(const Plus::procCMField<real> & parDiameter)override;

    /// Smooth vector field (no-op for PCM)
    void smoothenField(Foam::volVectorField& field)const override;
    
    /// Smooth scalar field (no-op for PCM)
    void smoothenField(Foam::volScalarField& field)const override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "PCM";
    }

};

}

#endif //__PCM_hpp__