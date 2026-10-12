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

#ifndef pFlow_multiReactiveThermalMomentumSphereUnresolvedCouplingSystem_hpp
#define pFlow_multiReactiveThermalMomentumSphereUnresolvedCouplingSystem_hpp

#include <vector>

#include "thermalMomentumSphereUnresolvedCouplingSystem.hpp"
#include "reactionInteraction.hpp"

namespace pFlow::coupling
{

/**
 * @brief Unresolved sphere coupling system with momentum, heat and
 * multi-species reaction coupling.
 *
 * The reaction sources of the DEM particles are mapped onto the CFD species
 * fields (speciesSu, speciesSp); the gas concentrations of the CFD cells are
 * sent back to the DEM particles.
 */
class multiReactiveThermalMomentumSphereUnresolvedCouplingSystem
:
    public thermalMomentumSphereUnresolvedCouplingSystem
{
private:

    /// Mapping of the reaction sources onto the mesh.
    reactionInteraction      reactionInteraction_;

    /// Explicit gas mass sources [kg/s], layout [p * N_gas + g].
    std::vector<real>        gasMassSource_;

    /// Picard implicit sink coefficients [kg.m3/(kmol.s)], same layout.
    std::vector<real>        gasMassSourceSp_;

    /// Gas concentrations [kmol/m3], layout [g * N_par + p], one segment per
    /// species.
    std::vector<real>        gasConcentrations_;

    /// Number of gas species.
    uint32                   numGasSpecies_ = 0;

    /// CFD species index of each DEM gas species.
    std::vector<int>         speciesIndexMap_;

    /// True until the DEM and CFD species have been cross-checked.
    bool                     firstMassCouplingCall_ = true;

    /// True when the DEM has gas species; set on all processors by
    /// buildSpeciesMapping(), as only the master processor holds the DEM.
    bool                     demHasGasSpecies_ = false;

    /// Gathers the gas mass sources and sinks of the DEM particles.
    bool collectReactionData();

    /// Sends the CFD gas concentrations to the DEM particles.
    bool dispatchGasConcentrationsToDEM();

    /**
     * Cross-checks the DEM and CFD gas species (names, order, molar masses)
     * once, in the first coupling step, and shares with all processors
     * whether the DEM has gas species. Stops with an error on a mismatch.
     */
    void buildSpeciesMapping();

public:

    TypeInfo(
        "multiReactiveSphereUnresolvedCouplingSystem"
        "<multiReactiveThermalMomentum>"
    );

    multiReactiveThermalMomentumSphereUnresolvedCouplingSystem
    (
        word            shapeTypeName,
        word            couplingSystemType,
        Foam::fvMesh&   mesh,
        int             argc,
        char* argv[]
    );

    multiReactiveThermalMomentumSphereUnresolvedCouplingSystem
    (
        const multiReactiveThermalMomentumSphereUnresolvedCouplingSystem&
    ) = delete;

    multiReactiveThermalMomentumSphereUnresolvedCouplingSystem& operator=
    (
        const multiReactiveThermalMomentumSphereUnresolvedCouplingSystem&
    ) = delete;

    ~multiReactiveThermalMomentumSphereUnresolvedCouplingSystem()
        override = default;

    add_vCtor
    (
        unresolvedCouplingSystem,
        multiReactiveThermalMomentumSphereUnresolvedCouplingSystem,
        word
    );

    void calculatePorosity() override;

    void calculateMassCoupling() override;

    bool sendDataToDEM(real t, real dt) override;

    /**
     * True when either the thermal/momentum base or the reaction
     * interaction needs the particle properties distributed with weights
     * over neighbouring cells instead of assigned to a single cell.
     */
    bool requireCellDistribution() const override
    {
        return
            thermalMomentumSphereUnresolvedCouplingSystem::
                requireCellDistribution()
         || reactionInteraction_.requireCellDistribution();
    }

    /// Explicit gas mass source of each species [kg/(m3.s)], for YEqn.
    inline const Foam::PtrList<Foam::volScalarField>& speciesSu() const override
    {
        return reactionInteraction_.speciesSu();
    }

    /// Implicit sink coefficient of each species [kg/(kmol.s)], for YEqn.
    inline const Foam::PtrList<Foam::volScalarField>& speciesSp() const override
    {
        return reactionInteraction_.speciesSp();
    }

    /// Explicit mass source [kg/(m3.s)] of the named gas species. Stops with
    /// an error if no species of that name is registered.
    Foam::tmp<Foam::volScalarField> massSource
    (
        const word& specieName
    ) const override;

    word couplingSystemType() const override
    {
        return "multiReactiveThermalMomentum";
    }
};

} // pFlow::coupling

#endif // pFlow_multiReactiveThermalMomentumSphereUnresolvedCouplingSystem_hpp
