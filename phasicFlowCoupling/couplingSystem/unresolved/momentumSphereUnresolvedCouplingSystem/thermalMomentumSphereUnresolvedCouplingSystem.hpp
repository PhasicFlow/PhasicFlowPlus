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

#ifndef pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp
#define pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp

#include "momentumSphereUnresolvedCouplingSystem.hpp"
#include "heatInteraction.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class thermalMomentumSphereUnresolvedCouplingSystem
 * @brief Adds heat transfer coupling on top of spherical momentum coupling.
 *
 * @details
 * Inherits porosity, momentum coupling (drag/lift/virtual mass) and their
 * source terms (Sp/Su/alpha) directly from
 * momentumSphereUnresolvedCouplingSystem - this class only adds
 * heatInteraction_ and the DEM<->CFD data exchange needed for heat
 * transfer (temperature, heat sources, and the PFP fluid property pipeline).
 */
class thermalMomentumSphereUnresolvedCouplingSystem
:
    public momentumSphereUnresolvedCouplingSystem
{
private:

    //- private members

        // --- Interaction models & managers ---

        /// @brief Convective + radiative heat interaction manager.
        heatInteraction             heatInteraction_;

        /// True when heatInteraction_ needs distribution weights beyond
        /// what the base class already provides.
        bool                        requiresDistributionHeat_ = false;

        // --- MPI communication fields (ProcCMFields) ---

        /// Particle temperature [K].
        Plus::realProcCMField       particleTemperature_;

        /// Convective heat source [W].
        Plus::realProcCMField       fluidHeatSourceConv_;

        /// Radiative heat source [W].
        Plus::realProcCMField       fluidHeatSourceRad_;

        /// Particle surface emissivity [-].
        Plus::realProcCMField       emissivity_;

        /// Sum of neighbouring particle temperatures [K].
        Plus::realProcCMField       radSumTemp_;

        /// Number of radiating neighbours [-].
        Plus::uint32ProcCMField     radNumPrt_;

        // PFP (Particle-Fluid-Particle) pipeline fields

        /// Local fluid thermal conductivity sampled at the particle [W/(m.K)].
        Plus::realProcCMField       fluidKappa_;

        /// Local fluid volume fraction sampled at the particle [-].
        Plus::realProcCMField       fluidAlpha_;

        // Diagnostic counters, accumulated over the whole run.

        /// Particle-timesteps with kappa floored to positive after
        /// MPI collection.
        uint64              cumulativeBadKappaCount_ = 0;

        /// Particle-timesteps with alpha clamped into [0,1] after
        /// MPI collection.
        uint64              cumulativeBadAlphaCount_ = 0;

        /// Particle-timesteps with an invalid mapped fluid cell index.
        uint64              cumulativeInvalidCellCount_ = 0;

        /// Next cumulative count at which a kappa/alpha reminder prints.
        uint64              nextKappaAlphaReportMilestone_ = 1;

        /// Next cumulative count at which an invalid-cell reminder prints.
        uint64              nextInvalidCellReportMilestone_ = 1;

    //- private methods

        // --- Private helper methods ---

        bool collectFluidHeatSource();
        
        bool collectFluidProperties();

        void sendFluidHeatSourceToDEM();
        
        void sendFluidPropertiesToDEM();

protected:

    //- protected methods

        /**
         * @brief Scatters DEM particle fields from the MPI master to 
         * all workers.
         */
        bool distributeParticleFields() override;

public:

    //- Type info

        TypeInfo(
            "thermalSphereUnresolvedCouplingSystem<thermalMomentum>");

    //- constructors

        // --- Constructors / destructor ---

        thermalMomentumSphereUnresolvedCouplingSystem(
            word            shapeTypeName,
            word            couplingSystemType,
            Foam::fvMesh&   mesh,
            int             argc,
            char*           argv[]);

        thermalMomentumSphereUnresolvedCouplingSystem(
            const thermalMomentumSphereUnresolvedCouplingSystem&) = delete;

        thermalMomentumSphereUnresolvedCouplingSystem& operator=(
            const thermalMomentumSphereUnresolvedCouplingSystem&) = delete;

        ~thermalMomentumSphereUnresolvedCouplingSystem() 
            override = default;

    //- public methods

        add_vCtor(
            unresolvedCouplingSystem,
            thermalMomentumSphereUnresolvedCouplingSystem,
            word
        );

        // --- Physical coupling calculations ---
        //
        // calculateMomentumCoupling(), calculateMassCoupling(), Sp(),
        // Su(), alpha() and shapeTypeName() are all inherited unchanged
        // from momentumSphereUnresolvedCouplingSystem: now that porosity_
        // and momentumInteraction_ live there, those methods already do
        // exactly what a thermal-specific override would duplicate.
        //
        // calculatePorosity() is the one exception - it still needs a
        // thin override (see the .C file): the base class only refreshes
        // distribution weights when ITS OWN requirement (porosity +
        // momentum) calls for it, which would silently skip the refresh
        // if only heatInteraction_ needs distribution.

        void calculatePorosity() override;
        
        void calculateHeatCoupling() override;

        // --- Source-term accessors ---

        Foam::tmp<Foam::volScalarField> heatSource() const override;

        /// @brief Implicit coefficient for the fluid energy equation 
        /// matrix [W/(m^3.K)].
        inline
        Foam::tmp<Foam::volScalarField> heatSp() const
        {
            return Foam::tmp<Foam::volScalarField>(
                heatInteraction_.heatSp());
        }

        /// @brief Explicit source term for the fluid energy equation 
        /// matrix [W/m^3].
        inline
        Foam::tmp<Foam::volScalarField> heatSu() const
        {
            return Foam::tmp<Foam::volScalarField>(
                heatInteraction_.heatSu());
        }

        // --- Identity ---

        inline
        word couplingSystemType() const override
        {
            return "thermalMomentum";
        }
        
        inline
        bool requireCellDistribution() const override
        {
            return 
                momentumSphereUnresolvedCouplingSystem::
                    requireCellDistribution() ||
                requiresDistributionHeat_;
        }

        // --- Data synchronisation ---

        bool sendDataToDEM(real t, real dt) override;

        inline
        Plus::realProcCMField& particleTemperature()
        {
            return particleTemperature_;
        }
        
        inline
        Plus::realProcCMField& fluidHeatSourceConv()
        {
            return fluidHeatSourceConv_;
        }
        
        inline
        Plus::realProcCMField& fluidHeatSourceRad()
        {
            return fluidHeatSourceRad_;
        }
        
        inline
        Plus::realProcCMField& fluidKappa()
        {
            return fluidKappa_;
        }
        
        inline
        Plus::realProcCMField& fluidAlpha()
        {
            return fluidAlpha_;
        }

        // --- Diagnostic accessors ---

        /// @brief Total particle-timesteps for which kappa or alpha needed
        /// clamping after MPI collection, accumulated over the whole run.
        inline
        uint64 cumulativeBadKappaCount() const
        {
            return cumulativeBadKappaCount_;
        }
        
        inline
        uint64 cumulativeBadAlphaCount() const
        {
            return cumulativeBadAlphaCount_;
        }

        /// @brief Total particle-timesteps for which the mapped fluid cell was
        /// invalid while sampling fluidKappa/fluidAlpha for the DEM side,
        /// accumulated over the whole run.
        inline
        uint64 cumulativeInvalidCellCount() const
        {
            return cumulativeInvalidCellCount_;
        }

}; // thermalMomentumSphereUnresolvedCouplingSystem

} // coupling
} // pFlow

#endif // pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp
