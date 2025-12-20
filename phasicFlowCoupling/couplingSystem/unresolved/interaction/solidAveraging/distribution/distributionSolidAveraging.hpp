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

/// @file distributionSolidAveraging.hpp
/// @brief Particle property averaging using spatial distribution kernels.
///
/// @class pFlow::coupling::distributionSolidAveraging
/// @ingroup interaction
/// @extends pFlow::coupling::solidAveraging
///
/// @brief Distribution-based property averaging with kernel smoothing.
///
/// Maps particle properties to fluid cells using distribution kernels
/// (Gaussian, diffusion, etc.), then back-interpolates to particles.
/// Ensures smooth spatial variation and reduces numerical noise.
///
/// ## Mathematical Formulation
///
/// **Step 1 - Volume Weighting:** Particle property weighted by volume
/// $$q_p^{vol} = V_p \cdot q_p$$
/// where $V_p = \frac{\pi}{6} d_p^3$ is particle volume.
///
/// **Step 2 - Distribution:** Weighted property distributed to cells
/// using kernel $W(d)$:
/// $$S_i = \sum_{p \in \text{kernel}(i)} W(d_{ip}) \cdot q_p^{vol}$$
///
/// **Step 3 - Normalization:** Divide by fluid volume in cell
/// $$\overline{q}_i = \frac{S_i}{(1-\alpha_i) V_i}$$
/// where $\alpha_i$ is solid volume fraction, $V_i$ is cell volume.
///
/// **Step 4 - Smoothing:** Apply filter (optional)
/// $$\overline{q}_i^{smooth} = \mathcal{F}(\overline{q}_i)$$
///
/// **Step 5 - Back-interpolation:** Return cell value at particle
/// $$\overline{q}_p = \overline{q}_{i(p)}$$
/// where $i(p)$ is cell containing particle $p$.
///
/// ## Data Members
/// - cellAvField_: Cell-centered averaged properties
/// - calcParAvField_: Back-interpolated particle values
///
/// @see solidAveraging (base class)
/// @see particleSolidAveraging (direct method)

#ifndef __distributionSolidAveraging_hpp__
#define __distributionSolidAveraging_hpp__

#include "solidAveraging.hpp"

namespace pFlow::coupling
{


class distributionSolidAveraging
:
    public solidAveraging
{
private:

    /// @brief Fluid field storing cell-centered averaged properties.
    Foam::volVectorField        cellAvField_;

    /// @brief Back-interpolated particle-centered averaged properties.
    Plus::realx3ProcCMField     calcParAvField_;

public:

    // type info
    TypeInfo("distribution");

    /// @brief Constructor for distribution-based property averaging.
    distributionSolidAveraging(
        const word&                     type,
        const unresolvedCouplingSystem& uCS,
        const porosity&                 prsty,
        const word&                     name
    );

    /// @brief Destructor.
    ~distributionSolidAveraging()override =default;

    /// @brief Virtual constructor macro for factory registration.
    add_vCtor
	(
		solidAveraging,
		distributionSolidAveraging,
		word	
	);
    
    /// @brief Perform mass-weighted averaging with distribution and smoothing.
    void calculate(
        const Plus::realx3ProcCMField& particleField) override;
    
    /// @brief Perform number-weighted averaging with distribution and smoothing.
    void calculateNumberBased(
        const Plus::realx3ProcCMField& particleField) override;
    
    /// @brief Cell distribution machinery is required for this method.
    bool requireCellDistribution()const override
    {
        return true;
    }

};

}

#endif //__distributionSolidAveraging_hpp__