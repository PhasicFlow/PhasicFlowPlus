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
 * @class GaussianIntegral
 * @brief Implements Gaussian integral method for void fraction calculation in particle-fluid coupling.
 * 
 * This class extends the distribution base class to provide Gaussian integral-based
 * distribution calculations for unresolved particle-fluid coupling systems. The method
 * uses a Gaussian kernel to distribute particle properties to the computational mesh.
 * 
 * @details The implementation is based on the research by Alireza Kianimoqadam and Justin Lapp,
 * "Gaussian integral method for void fraction", Particuology, Volume 108, 2026, Pages 125-142.
 * 
 * The class maintains neighbor lists within a specified radius (maxLayers_) around each
 * cell to efficiently compute the Gaussian integral distribution. Lists are constructed
 * lazily and cached for performance.
 * 
 * @see distribution
 * @see couplingMesh
 * @see Plus::centerMassField
 */

 
#ifndef __GaussianIntegral_hpp__ 
#define __GaussianIntegral_hpp__


#include "distribution.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class GaussianIntegral
:
    public distribution
{
private:
    
    // radius of circule for cell neighbor list 
    Foam::label 				maxLayers_;

    bool                		listsConstructed_ = false;

public:

    /// Type info
    TypeInfo("GaussianIntegral");

    /// Construc from dictionary 
    GaussianIntegral(
        Foam::dictionary 		dict, 
        const couplingMesh& 	cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    ~GaussianIntegral() = default;

    add_vCtor
    (
        distributionBase,
        GaussianIntegral,
        dictionary  
    );

    void checkForListsConstructed();

    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    Foam::word distributionMethodName()const override
    {
        return "GaussianIntegral";
    }

}; 

} // pFlow::coupling


#endif
