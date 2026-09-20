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

#ifndef pFlow_sphereHeatTransfer_hpp 
#define pFlow_sphereHeatTransfer_hpp

#include "heatTransfer.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "distributionBase.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class sphereHeatTransfer
 * @brief Computes convective and radiative heat transfer for spherical 
 *        particles.
 *
 * @details 
 * This class executes the core physics of heat exchange between the 
 * fluid mesh and the discrete spherical particles. 
 * - Convection: Calculated via the injected `NusseltClosureType` 
 *   (e.g., Ranz-Marshall).
 *   The resulting source terms (Su, Sp) are the ONLY contribution this class
 *   makes to the fluid energy equation.
 * - Radiation: Implements a neighbourhood-averaged radiation model between
 *   solid particles. The carrier fluid is treated as radiatively transparent
 *   - it never emits or absorbs radiation, so radiative exchange affects
 *   only QpRad (the heat absorbed by each particle) and is deliberately kept
 *   out of Su/Sp. The "environment" temperature driving this exchange blends
 *   the local fluid temperature and the neighbourhood-averaged particle
 *   temperature purely as a numerical proxy for particles that are not
 *   individually resolved within the radiation cut-off radius; it does not
 *   imply the fluid itself participates energetically.
 *
 * @note This class uses C++ Templates to resolve the Nusselt model at 
 * compile-time, avoiding the overhead of virtual function calls inside the 
 * massive particle loop.
 */
template<typename NusseltClosureType>
class sphereHeatTransfer 
: 
    public heatTransfer
{
public:

    //- Type info

        /// @brief Alias for the template type to simplify macro registration.
        using SphereHeatTransferType = sphereHeatTransfer<NusseltClosureType>;

        TypeInfoTemplate11("sphereHeatTransfer", NusseltClosureType);

private:

    //- private members

        // --- Closure models and state flags ---

        /// @brief The specific closure model (functor) for calculating the 
        /// Nusselt number.
        NusseltClosureType  nusseltClosure_;

        /// @brief Flag to ensure the unresolved assumption warning is printed 
        /// only once.
        bool                unresolvedWarningIssued_ = false;

public:

    //- constructors

        // --- Constructors & destructor ---

        sphereHeatTransfer(
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty);

        virtual ~sphereHeatTransfer() override = default;

        add_vCtor(
            heatTransfer, 
            SphereHeatTransferType, 
            couplingSystem
        );

    //- public methods

        // --- Core physics routine ---

        /**
         * @brief Core mathematical function to calculate convection and 
         * radiation forces.
         * @details Evaluates particle heat fluxes and distributes ONLY the
         * convective source terms (Su, Sp) back onto the Eulerian fluid mesh.
         * Radiative exchange is computed as well, but is returned exclusively
         * through QpRad and never contributes to Su/Sp - the fluid is
         * radiatively transparent and must have zero net energy exchange
         * through this mechanism. See the "Radiation Calculation" comment
         * block inside calculateHeatTransfer(), in sphereHeatTransfer.C,
         * for the physical modelling assumptions behind this split.
         * 
         * @param fluidVelocity Interpolated fluid velocity at particle centers.
         * @param parVelocity Discrete particle velocity.
         * @param diameter Discrete particle diameter.
         * @param parTemperature Discrete particle temperature.
         * @param cellDistribution Operator mapping discrete data to mesh cells.
         * @param Qp [OUT] Convective heat flux injected into particles.
         * @param emissivity Discrete particle emissivity.
         * @param radSumTemp Sum of neighboring particle temperatures.
         * @param radNumPrt Number of neighbors contributing to radiation.
         * @param QpRad [OUT] Radiative heat flux injected into particles. This
         * is the ONLY output of the radiation calculation; it is consumed by
         * the DEM-side particle temperature integrator and never reaches the
         * fluid energy equation.
         */
        void calculateHeatTransfer(
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    diameter,
            const Plus::realProcCMField&    parTemperature,
            const distributionBase&         cellDistribution,  
            Plus::realProcCMField&          Qp,                    
            const Plus::realProcCMField&    emissivity,            
            const Plus::realProcCMField&    radSumTemp,            
            const Plus::uint32ProcCMField&  radNumPrt,             
            Plus::realProcCMField&          QpRad) override;

}; 

} // coupling
} // pFlow

// Include the template implementation file
#include "sphereHeatTransfer.C"

#endif // pFlow_sphereHeatTransfer_hpp
