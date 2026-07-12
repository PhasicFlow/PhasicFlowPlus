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
 * @class adaptiveGaussian
 * @brief Adaptive Gaussian distribution method for particle-to-cell coupling.
 *
 * @details
 * This class implements an adaptive Gaussian distribution method that maps particle
 * properties to fluid cells using a Gaussian kernel with adaptive standard deviation.
 * The standard deviation is dynamically adjusted based on the particle-to-cell size 
 * ratio, making this method suitable for unresolved Eulerian-Lagrangian multiphase 
 * flow simulations.
 *
 * **Mathematical Formulation:**
 *
 * For a particle at position $\mathbf{x}_p$ with property value $q_p$, the weight 
 * assigned to cell $i$ is calculated using a Gaussian kernel:
 *
 * $$w_i = \frac{\exp\left(-\frac{|\mathbf{x}_p - \mathbf{x}_i|^2}{2\sigma^2}\right)}{\sum_j \exp\left(-\frac{|\mathbf{x}_p - \mathbf{x}_j|^2}{2\sigma^2}\right)}$$
 *
 * The adaptive standard deviation is computed as:
 *
 * $$\sigma = d_{cell} \cdot f_s \cdot a \cdot \left(\frac{d_{cell}}{d_p}\right)^{e}$$
 *
 * where:
 * - $d_{cell} = V_{cell}^{1/3}$ is the characteristic cell size
 * - $f_s$ is a smoothing factor (user-configurable, default: 1.0)
 * - $a = 0.6142275$ is an empirical constant
 * - $e = -0.6195039$ is an empirical exponent
 * - $d_p$ is the particle diameter
 *
 * **Key Features:**
 *
 * - **Automatic mesh adaptation**: The standard deviation adapts to different cell 
 *   sizes and particle diameters, making the method robust across various mesh resolutions.
 * - **High accuracy**: The method is accurate for cell-to-particle size ratios down to 1, 
 *   making it suitable for nearly all unresolved simulations.
 * - **Computational efficiency**: Compared to sub-division methods, it provides good 
 *   computational performance.
 * - **Smooth distributions**: Produces smooth transitions of particle properties across 
 *   neighboring cells.
 * - **Efficient neighbor search**: Uses neighbor lists with configurable maximum layers 
 *   to limit the search radius and improve computational efficiency.
 *
 * **Dictionary Configuration:**
 *
 * ```cpp
 * adaptiveGaussianInfo
 * {
 *     maxLayers           1;          // optional, default: 1
 *     smoothingFactor     1.0;        // optional, default: 1.0
 * }
 * ```
 *
 * @param maxLayers Maximum number of cell layers to consider for neighbor search 
 *                  (default: 1)
 * @param smoothingFactor User-defined scaling factor for standard deviation 
 *                        (default: 1.0)
 *
 * @see distribution
 * @see distributionBase
 *
 * @note 
 * For cell-to-particle size ratios greater than 7, the method automatically reverts 
 * to Particle Centroid Method (PCM) for computational efficiency.
 */

#ifndef __adaptiveGaussian_hpp__
#define __adaptiveGaussian_hpp__

#include "distribution.hpp"

namespace pFlow::coupling
{

class couplingMesh;


class adaptiveGaussian
:
    public distribution
{
    /// Maximum number of cell layers for neighbor search
    Foam::label         maxLayers_;

    /// Flag for neighbor list construction state
    bool                listsConstructed_ = false;

    /// Empirical constant for standard deviation calculation
    const static inline Foam::scalar   a_ = 0.6142275;

    /// Empirical exponent for standard deviation calculation
    const static inline Foam::scalar   exponent_ = -0.6195039;

    /// User-defined smoothing factor for standard deviation
    static inline Foam::scalar         smoothingFactor_ = 1.0;

    /// Calculate adaptive standard deviation based on particle-to-cell size ratio
    inline 
    Foam::scalar stdDeviation(const Foam::scalar& dx_dp, const Foam::scalar& dcell) const
    {
        return dcell* smoothingFactor_* a_ * Foam::pow(dx_dp, exponent_);
    }

    

public:

    /// Type info
    TypeInfo("adaptiveGaussian");

    /// Initialize the adaptive Gaussian distribution from dictionary
    adaptiveGaussian(
        const Foam::dictionary& 	 parrentDict, 
        const couplingMesh& 	     cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    virtual ~adaptiveGaussian() = default;

    add_vCtor
    (
        distributionBase,
        adaptiveGaussian,
        dictionary
    );

    /// Construct neighbor lists lazily on first call
    void checkForListsConstructed();

    /// Update distribution weights based on particle diameters
    void updateWeights(const Plus::procCMField<real> & parDiameter) override;

    /// Return the name of the distribution method
    Foam::word distributionMethodName()const override
    {
        return "adaptiveGaussian";
    }

    

}; // end adaptiveGaussian

} // end pFlow::coupling


#endif //__adaptiveGaussian_hpp__
