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
 * @file diffusion.hpp
 * @class diffusion
 * @brief Laplacian diffusion distribution method for particle-to-cell coupling.
 *
 * @details
 * Implements a diffusion-based distribution method that smooths particle properties 
 * across fluid cells using the Laplacian diffusion equation. Starting from Particle 
 * Centroid Method (PCM) initial conditions, the method applies iterative diffusion 
 * steps to produce smooth field distributions suitable for unresolved multiphase 
 * flow simulations.
 *
 * **Mathematical Formulation:**
 *
 * The diffusion equation is solved iteratively:
 *
 * $$\frac{\partial \phi}{\partial \tau} = D \nabla^2 \phi$$
 *
 * where:
 * - $D = \sigma^2 / (4 \tau_{total})$ is the diffusion coefficient
 * - $\sigma$ is the standard deviation parameter (user-configurable)
 * - $\tau_{total} = 1 \text{ s}$ is the total pseudo-time
 * - Time step: $\Delta \tau = \tau_{total} / nSteps$
 *
 * **Key Features:**
 *
 * - **Smooth field transitions:** Produces physically smoother distributions across neighboring cells
 * - **Tunable smoothing:** Control over smoothness via `standardDeviation` and `nSteps` parameters
 * - **Time-convergent:** Results are reproducible and converge to steady-state distribution
 * - **Good accuracy-cost balance:** Computationally efficient relative to sub-division methods
 * - **Flexible configuration:** Optional solver output to screen via `log` parameter
 *
 * **Dictionary Configuration:**
 *
 * ```cpp
 * diffusionInfo
 * {
 *     nSteps              5;              // Number of diffusion steps (required)
 *     standardDeviation   0.0075;         // Standard deviation parameter (required)
 *     log                 0;              // Optional, log solver output (default: 0)
 * }
 * ```
 *
 * @param nSteps Number of pseudo-time diffusion steps to apply
 * @param standardDeviation Standard deviation controlling diffusion extent
 * @param log Optional flag to enable solver output (default: 0)
 *
 * @see distributionBase
 * @see adaptiveGaussian
 *
 * @note
 * Recommended `standardDeviation` values are 3-6 times the particle diameter, 
 * depending on the cell-to-particle size ratio.
 */

#ifndef __diffusion_hpp__
#define __diffusion_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"


// from phasicFlowPlus
#include "distributionBase.hpp"


namespace pFlow::coupling
{


class diffusion
:
    public distributionBase
{
    /// Number of pseudo-time diffusion steps
    Foam::label 			    nSteps_;

	/// Standard deviation controlling diffusion extent
	Foam::scalar 			    standardDeviation_;

	/// Total integration time for diffusion
	Foam::scalar 			    intTime_;

	/// Time step for each diffusion iteration
	Foam::dimensionedScalar     dt_;

	/// Diffusion coefficient
	Foam::dimensionedScalar     DT_;

	/// Dictionary for smooth field solver configuration
	Foam::dictionary            smoothSolDict_;

    /// Construct time derivative matrix for scalar fields
    Foam::tmp<Foam::fvMatrix<Foam::scalar>> fvmDdt
    (
        const Foam::volScalarField& sField
    )const;

    /// Construct time derivative matrix for vector fields
    Foam::tmp<Foam::fvMatrix<Foam::vector>> fvmDdt
    (
        const Foam::volVectorField& sField
    )const;

protected:
    
    /// Construct neighbor lists with specified search length and max layers
    void constructLists(
        const Foam::scalar searchLen, 
        const Foam::label maxLayers) override;
    
    /// Construct neighbor lists with max layers only
    void constructLists(const Foam::label maxLayers)override;

public:

    /// Type info
	TypeInfo("diffusion");

	/// Initialize diffusion distribution from dictionary
	diffusion(
		const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~diffusion()=default;
    
    add_vCtor
    (
        distributionBase,
        diffusion,
        dictionary  
    );

    /// Update distribution weights through iterative diffusion steps
    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    /// Apply diffusion smoothing to scalar field
    void smoothenField(Foam::volScalarField& field)const override;

    /// Apply diffusion smoothing to vector field
    void smoothenField(Foam::volVectorField& field)const override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "diffusion";
    }

};

}

#endif //__diffusion_hpp__
