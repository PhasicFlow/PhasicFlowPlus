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

#ifndef pFlow_heatInteraction_hpp
#define pFlow_heatInteraction_hpp

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "Timer.hpp"
#include "procCMFields.hpp"
#include "heatTransfer.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "PCM.hpp"

namespace pFlow
{
namespace coupling
{

// Forward declarations
class unresolvedCouplingSystem;
class porosity;

/**
 * @class heatInteraction
 * @brief Manages fluid-particle heat transfer, mirroring momentumInteraction.
 *
 * @details
 * `fluidVelocity`/`solidVelocity` may each independently be set to
 * "similarToMomentum" to reuse momentumInteraction's averaging instead of
 * building a dedicated one for heat. `heatSourceExchange` separately
 * selects how the heat source is mapped onto the fluid mesh
 * ("distribution" or "cell") - it has no "similarToMomentum" option.
 */
class heatInteraction
{
private:

    //- private members

        /// Reference to the porosity field for cells
        const porosity&             porosity_;

        /// Heat-only fluid velocity averaging; null when "fluidVelocity"
        /// is "similarToMomentum" (calculateCoupling() then falls back
        /// to the reference it is called with)
        uniquePtr<fluidAveraging>   fluidAveraging_ = nullptr;

        /// Heat-only solid velocity averaging; null when "solidVelocity"
        /// is "similarToMomentum" (calculateCoupling() then falls back
        /// to the reference it is called with)
        uniquePtr<solidAveraging>   solidAveraging_ = nullptr;

        /// Flag indicating if the heat source is exchanged using
        /// distribution (true) or cell method (false)
        bool                        heatExchangeDistribute_;
        
        /// Flag indicating if cell distribution is required for the coupling
        bool                        requireCellDistribution_ = false;

        /// Pointer to the heat transfer closure model (e.g. Ranz-Marshall)
        uniquePtr<heatTransfer>     heatTransfer_;
        
        /// Pointer to the particle-centroid mapper for non-distributed
        /// heat exchange
        uniquePtr<PCM>              noDistribution_ = nullptr;
        
        /// Timer for tracking performance of heat interaction calculations
        Timer                       heatInteractionTimer_;

public:
    
    //- Type info

        TypeInfo("heatInteraction");

    //- constructors

        heatInteraction(
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty);
        
        virtual ~heatInteraction() = default;

    //- public methods

        /// @brief Implicit source coefficient (Sp) for the fluid energy 
        /// equation matrix.
        inline
        const Foam::volScalarField& heatSp() const 
        { 
            return std::as_const<const heatTransfer&>(*heatTransfer_).Sp(); 
        }

        /// @brief Explicit source vector (Su) for the fluid energy equation 
        /// matrix.
        inline
        const Foam::volScalarField& heatSu() const 
        { 
            return std::as_const<const heatTransfer&>(*heatTransfer_).Su(); 
        }

        /// @brief Alias for heatSu, returning the total volumetric heat source.
        inline
        const Foam::volScalarField& Sh() const 
        { 
            return heatSu(); 
        }

        inline
        const porosity& Porosity() const
        {
            return porosity_;
        }

        const unresolvedCouplingSystem& uCS() const;

        inline
        bool requireCellDistribution() const 
        { 
            return requireCellDistribution_; 
        }

        const Foam::dictionary& dict() const;

        static const Foam::dictionary& getDict(
            const unresolvedCouplingSystem& uCS);

        /**
         * @brief Calculates convective and radiative heat exchanges.
         * @details Gathers required physical fields and delegates the actual 
         * calculation to the instantiated `heatTransfer` model, applying the 
         * user-selected distribution mapping. If this object owns its own
         * fluidAveraging_/solidAveraging_ (i.e. the dictionary did not say
         * "similarToMomentum" for that key), they are calculated here from
         * U/vp and used in place of the fluidVelocity/parVelocity arguments;
         * otherwise fluidVelocity/parVelocity (momentum's already-computed
         * averaging) are used directly and U/vp go unused.
         *
         * @param U             Fluid velocity field - only used when
         *                      fluidAveraging_ is non-null.
         * @param vp            Particle velocity - only used when
         *                      solidAveraging_ is non-null.
         * @param fluidVelocity Momentum's fluid velocity averaging - used
         *                      only when fluidAveraging_ is null.
         * @param parVelocity   Momentum's solid velocity averaging - used
         *                      only when solidAveraging_ is null.
         * @param dp            Particle diameter.
         * @param Tp            Particle temperature.
         * @param Qp            [OUT] Computed convective heat source.
         * @param emissivity    Particle emissivity.
         * @param radSumTemp    Sum of neighboring temperatures for radiation.
         * @param radNumPrt     Number of neighbors considered in radiation.
         * @param QpRad         [OUT] Computed radiative heat source.
         */
        virtual void calculateCoupling(
            const Foam::volVectorField&     U,
            const Plus::realx3ProcCMField&  vp,
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    dp,
            const Plus::realProcCMField&    Tp,
            Plus::realProcCMField&          Qp,
            const Plus::realProcCMField&    emissivity,
            const Plus::realProcCMField&    radSumTemp,
            const Plus::uint32ProcCMField&  radNumPrt,
            Plus::realProcCMField&          QpRad);

}; // heatInteraction

} // coupling
} // pFlow

#endif // pFlow_heatInteraction_hpp
