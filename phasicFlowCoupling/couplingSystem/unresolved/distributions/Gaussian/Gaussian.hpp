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
 * @class Gaussian
 * @brief Gaussian distribution method for particle-to-cell coupling.
 *
 * Distributes particle properties to fluid cells using a Gaussian kernel 
 * with user-specified standard deviation and configurable neighbor search radius.
 */

#ifndef __Gaussian_hpp__ 
#define __Gaussian_hpp__

// from PhasicFlowPlus
#include "distribution.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class Gaussian
:
    public distribution
{	
private:
    
    /// Standard deviation of the Gaussian kernel
    Foam::scalar 				standardDeviation_;

    /// Maximum number of cell layers for neighbor search
    Foam::label 				maxLayers_;

    /// List of boundary cell pairs
    std::vector<std::pair<Foam::label, Foam::label>>    boundaryCell_;

    /// Construct boundary cell lists with specified search length
    void constructBoundaryLists(const Foam::scalar searchLen);

    /// Flag for neighbor list construction state
    bool                		listsConstructed_ = false;

public:

    /// Type info
    TypeInfo("Gaussian");

    /// Initialize Gaussian distribution from dictionary
    Gaussian(
        Foam::dictionary 		dict, 
        const couplingMesh& 	cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    ~Gaussian() = default;

    add_vCtor
    (
        distributionBase,
        Gaussian,
        dictionary  
    );

    /// Get the standard deviation of the Gaussian kernel
    inline
    auto standardDeviation()const
    {
        return standardDeviation_;
    }

    /// Construct neighbor lists lazily on first call
    void checkForListsConstructed();

    /// Update distribution weights based on particle diameters
    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "Gaussian";
    }
}; 

} // pFlow::coupling


#endif
