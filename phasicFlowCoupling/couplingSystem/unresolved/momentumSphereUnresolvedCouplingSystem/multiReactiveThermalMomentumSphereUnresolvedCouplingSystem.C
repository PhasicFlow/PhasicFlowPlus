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

#include <algorithm>
#include <string>

#include "multiReactiveThermalMomentumSphereUnresolvedCouplingSystem.hpp"
#include "procCommunicationPlus.hpp"

namespace pFlow::coupling
{

//------------------------------ private methods -------------------------------

bool multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
collectReactionData()
{
    auto& comm = parMapping().realScatteredComm();
    const size_t N_gas  = static_cast<size_t>(numGasSpecies_);
    const size_t nLocal = this->parMapping().centerMass().size();

    // If CFD has no species, or the DEM has multi-reaction disabled, skip
    // the MPI distribution loop. This allows pure-thermal runs.
    // demHasGasSpecies_ is the same on every processor (see
    // buildSpeciesMapping()), so all of them take the same path.
    if (N_gas == 0 || !demHasGasSpecies_) return true;

    auto allGasSrc = pDEMSystem().gasMassSource();
    auto allGasSp  = pDEMSystem().gasMassSourceSp();

    const size_t nTotal =
        (allGasSrc.size() > 0) ? allGasSrc.size() / N_gas : size_t(0);

    std::vector<real> globalSpec
    (
        Plus::processor::isMaster() ? nTotal : size_t(0),
        real(0)
    );
    std::vector<real> localSpec(nLocal, real(0));

    // Scatters one species column at a time, so the same routine handles any
    // number of gas species.
    auto scatterSpeciesData =
        [&](const auto& allData, std::vector<real>& localDst) -> bool
    {
        for (size_t g = 0; g < N_gas; ++g)
        {
            if (Plus::processor::isMaster() && nTotal > 0)
            {
                for (size_t p = 0; p < nTotal; ++p)
                {
                    globalSpec[p] = allData[p * N_gas + g];
                }
            }

            pFlow::span<real> globalSpan(globalSpec.data(), globalSpec.size());
            pFlow::span<real> localSpan (localSpec.data(), localSpec.size());

            if (!comm.distribute(globalSpan, localSpan)) return false;

            for (size_t p = 0; p < nLocal; ++p)
            {
                localDst[p * N_gas + g] = localSpec[p];
            }
        }
        return true;
    };

    if (!scatterSpeciesData(allGasSrc, gasMassSource_)) return false;
    if (!scatterSpeciesData(allGasSp, gasMassSourceSp_)) return false;

    return true;
}

bool multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
dispatchGasConcentrationsToDEM()
{
    // No CFD gas species: nothing to send. numGasSpecies_ is the same on
    // every processor.
    if (numGasSpecies_ == 0) return true;

    if (firstMassCouplingCall_)
    {
        buildSpeciesMapping();
        firstMassCouplingCall_ = false;
    }

    // No DEM species (reactions disabled): nothing to send
    if (!demHasGasSpecies_) return true;

    const auto& cellIDs  = this->parCellIndex();
    const size_t numPar  = cellIDs.size();
    const size_t N_gas   = static_cast<size_t>(numGasSpecies_);
    const size_t reqSize = numPar * N_gas;

    if (gasConcentrations_.size() != reqSize)
    {
        gasConcentrations_.assign(reqSize, real(0));
    }
    else
    {
        std::fill
        (
            gasConcentrations_.begin(),
            gasConcentrations_.end(),
            real(0)
        );
    }

    const auto& rho =
        this->cMesh().mesh().lookupObject<Foam::volScalarField>("rho");
    const auto& tp =
        this->cMesh().mesh().lookupObject<Foam::IOdictionary>
        ("transportProperties");
    Foam::List<Foam::scalar> gasMw(tp.lookup("gasMw"));

    // Concentration of species g in the cell of particle p: rho*Y/Mw
    for (size_t g = 0; g < N_gas; ++g)
    {
        const int cfdIdx = speciesIndexMap_[g];
        const auto& Y_g =
            this->cMesh().mesh().lookupObject<Foam::volScalarField>
            ("Y." + reactionInteraction_.gasSpeciesNames()[g]);

        for (size_t p = 0; p < numPar; ++p)
        {
            const Foam::label cellI = cellIDs[p];
            gasConcentrations_[g * numPar + p] =
                (cellI >= 0 && cellI < static_cast<Foam::label>(rho.size()))
                ? static_cast<real>(rho[cellI] * Y_g[cellI])
                    / static_cast<real>(gasMw[cfdIdx])
                : real(0);
        }
    }

    auto& comm        = parMapping().realScatteredComm();
    auto  allGasConcs = pDEMSystem().gasConcentrations();
    const size_t nTotal =
        Plus::processor::isMaster() && allGasConcs.size() > 0
            ? allGasConcs.size() / N_gas
            : size_t(0);

    if (Plus::processor::isMaster() && nTotal > 0)
    {
        std::fill
        (
            allGasConcs.data(),
            allGasConcs.data() + allGasConcs.size(),
            real(0)
        );
    }

    std::vector<real> globalSpec
    (
        Plus::processor::isMaster() ? nTotal : size_t(0),
        real(0)
    );

    for (size_t g = 0; g < N_gas; ++g)
    {
        // collectSum() adds into the receive buffer, so clear it for each
        // species
        std::fill(globalSpec.begin(), globalSpec.end(), real(0));

        // Each species is sent from its own segment of the member buffer:
        // collectSum() does not wait for the sends to complete
        pFlow::span<real> localSpan
        (
            gasConcentrations_.data() + g * numPar,
            numPar
        );
        pFlow::span<real> globalSpan(globalSpec.data(), globalSpec.size());

        if (!comm.collectSum(localSpan, globalSpan)) return false;

        if (Plus::processor::isMaster() && nTotal > 0)
        {
            for (size_t p = 0; p < nTotal; ++p)
            {
                allGasConcs[p * N_gas + g] = globalSpec[p];
            }
        }
    }

    return pDEMSystem().sendGasConcentrationsToDEM();
}

void multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
buildSpeciesMapping()
{
    const std::vector<std::string> stdNames = pDEMSystem().gasSpeciesNames();

    {
        Plus::procCommunication proc;

        real hasSpeciesMaster =
            (Plus::processor::isMaster() && stdNames.size() > 0)
          ? real(1)
          : real(0);
        real hasSpecies = real(0);

        if (!proc.distributeMasterToAll(hasSpeciesMaster, hasSpecies))
        {
            fatalErrorInFunction
                << "could not share the DEM gas species flag.\n";
            Plus::processor::abort(0);
        }

        demHasGasSpecies_ = hasSpecies > real(0.5);
    }

    if (!demHasGasSpecies_ && numGasSpecies_ > 0)
    {
        Foam::Info
            << "\n[CLC] The DEM has no gas species (multi-reaction disabled):"
            << " reaction sources stay zero and no gas concentrations are "
            << "sent to the DEM.\n" << Foam::endl;
    }

    const auto& msh = this->cMesh().mesh();

    if (!msh.foundObject<Foam::IOdictionary>("transportProperties"))
    {
        FatalErrorInFunction
            << "transportProperties not found in the OpenFOAM ObjectRegistry."
            << Foam::exit(Foam::FatalError);
    }

    const auto& tp =
        msh.lookupObject<Foam::IOdictionary>("transportProperties");
    Foam::wordList tpSpecies(tp.lookup("gasSpecies"));
    Foam::List<Foam::scalar> tpGasMw(tp.lookup("gasMw"));

    if (tpGasMw.size() != tpSpecies.size())
    {
        FatalErrorInFunction
            << "transportProperties: 'gasMw' has " << tpGasMw.size()
            << " entries but 'gasSpecies' has " << tpSpecies.size()
            << " entries. The two lists must be the same length and in "
            << "the same order." << Foam::exit(Foam::FatalError);
    }

    speciesIndexMap_.resize(numGasSpecies_, -1);

    Foam::Info << "\n[CLC] Building species index map:" << Foam::endl;

    for (size_t g = 0; g < numGasSpecies_; ++g)
    {
        const Foam::word& cpName = reactionInteraction_.gasSpeciesNames()[g];

        bool found = false;
        for (Foam::label i = 0; i < tpSpecies.size(); ++i)
        {
            if (cpName == tpSpecies[i])
            {
                speciesIndexMap_[g] = i;
                found = true;
                Foam::Info
                    << "  [" << g << "] " << cpName << " -> " << i
                    << Foam::endl;
                break;
            }
        }

        if (!found)
        {
            FatalErrorInFunction
                << "Species '" << cpName << "' NOT FOUND."
                << Foam::exit(Foam::FatalError);
        }
    }

    // The DEM comparisons below use master-only data
    if (!Plus::processor::isMaster() || !demHasGasSpecies_) return;

    Foam::wordList demNames(static_cast<Foam::label>(stdNames.size()));
    for (size_t i = 0; i < stdNames.size(); ++i)
    {
        demNames[i] = Foam::word(stdNames[i]);
    }

    if (static_cast<size_t>(demNames.size()) != numGasSpecies_)
    {
        FatalErrorInFunction
            << "Species count mismatch between DEM and CFD coupling!"
            << Foam::exit(Foam::FatalError);
    }

    // Identity check: DEM and CFD must name the same gas species in the same
    // order used for gasConcentrations() and gasMassSource(). speciesIndexMap_
    // ties DEM index g to a CFD species by name; this also verifies that the
    // DEM list, demNames[g], names that same species at that position.
    for (size_t g = 0; g < numGasSpecies_; ++g)
    {
        const Foam::word demSpecies(demNames[g]);
        const Foam::word cfdSpecies = tpSpecies[speciesIndexMap_[g]];

        if (demSpecies != cfdSpecies)
        {
            FatalErrorInFunction
                << "Gas species order mismatch between DEM and CFD at index "
                << g
                << ": DEM reports '" << demSpecies
                << "', CFD (transportProperties/gasSpecies) maps to '"
                << cfdSpecies << "'.\n"
                << "The species and their order must agree across "
                << "transportProperties (gasSpecies), couplingProperties "
                << "(reactionInteraction/reactionTransferModel/productGas), "
                << "and thermoPhysicalInteraction "
                << "(multiReactiveProperties/productGas)."
                << Foam::exit(Foam::FatalError);
        }
    }

    // Molar mass cross-check: catches a transportProperties/gasMw entry that
    // differs from the DEM value of a species whose name and order agree.
    const std::vector<real> demGasMw = pDEMSystem().gasMolarMasses();

    if (demGasMw.size() > 0)
    {
        if (demGasMw.size() != numGasSpecies_)
        {
            FatalErrorInFunction
                << "DEM reports " << demGasMw.size()
                << " gas molar masses but " << numGasSpecies_
                << " gas species names." << Foam::exit(Foam::FatalError);
        }

        constexpr Foam::scalar relTol = 1e-3;

        for (size_t g = 0; g < numGasSpecies_; ++g)
        {
            const Foam::scalar demMw = static_cast<Foam::scalar>(demGasMw[g]);
            const Foam::scalar cfdMw = tpGasMw[speciesIndexMap_[g]];

            const Foam::scalar denom = Foam::max(Foam::mag(demMw), Foam::SMALL);
            const Foam::scalar relDiff = Foam::mag(demMw - cfdMw) / denom;

            if (relDiff > relTol)
            {
                const Foam::word speciesName(demNames[g]);
                FatalErrorInFunction
                    << "Gas molar mass mismatch for species '"
                    << speciesName << "': "
                    << "DEM (thermoPhysicalInteraction/gasMolarMasses) = "
                    << demMw << " kg/kmol, CFD (transportProperties/gasMw) = "
                    << cfdMw << " kg/kmol (relative difference "
                    << relDiff * 100.0 << "%, "
                    << "tolerance " << relTol * 100.0 << "%).\n"
                    << "Update whichever dictionary is stale so both sides "
                    << "agree."
                    << Foam::exit(Foam::FatalError);
            }
        }
    }
}

//-------------------------------- constructors --------------------------------

multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
multiReactiveThermalMomentumSphereUnresolvedCouplingSystem
(
    word            shapeTypeName,
    word            couplingSystemType,
    Foam::fvMesh&   mesh,
    int             argc,
    char* argv[]
)
:
    thermalMomentumSphereUnresolvedCouplingSystem
    (
        shapeTypeName,
        couplingSystemType,
        mesh,
        argc,
        argv
    ),
    reactionInteraction_(*this, this->porosityCoupling()),
    gasMassSource_  (),
    gasMassSourceSp_()
{
    numGasSpecies_ =
        static_cast<uint32>(reactionInteraction_.speciesSu().size());

    Foam::Info
        << "\n[CLC Extension] Multi-Reactive Coupling System initialised.\n"
        << "    -> Mapped gas species count: " << numGasSpecies_
        << "\n" << Foam::endl;
}

//------------------------------- public methods -------------------------------

void multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
calculatePorosity()
{
    thermalMomentumSphereUnresolvedCouplingSystem::calculatePorosity();

    // Refreshes the distribution weights when only the reaction interaction
    // needs them
    if
    (
        reactionInteraction_.requireCellDistribution()
     && !thermalMomentumSphereUnresolvedCouplingSystem::
            requireCellDistribution()
    )
    {
        this->updateDistributionWeights();
    }
}

void multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
calculateMassCoupling()
{
    if (firstMassCouplingCall_)
    {
        buildSpeciesMapping();
        firstMassCouplingCall_ = false;
    }

    const size_t nLocal       = this->parMapping().centerMass().size();
    const size_t requiredSize = nLocal * static_cast<size_t>(numGasSpecies_);

    if (gasMassSource_.size() != requiredSize)
    {
        gasMassSource_.assign(requiredSize, real(0));
    }
    else
    {
        std::fill(gasMassSource_.begin(), gasMassSource_.end(), real(0));
    }

    if (gasMassSourceSp_.size() != requiredSize)
    {
        gasMassSourceSp_.assign(requiredSize, real(0));
    }
    else
    {
        std::fill(gasMassSourceSp_.begin(), gasMassSourceSp_.end(), real(0));
    }

    if (!collectReactionData())
    {
        fatalErrorInFunction << "collectReactionData() failed.\n";
        Plus::processor::abort(0);
    }

    reactionInteraction_.calculateCoupling(gasMassSource_, gasMassSourceSp_);
}

bool multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::
sendDataToDEM(real t, real dt)
{
    if (!thermalMomentumSphereUnresolvedCouplingSystem::sendDataToDEM(t, dt))
    {
        return false;
    }

    sendDataTimer().start();

    if (!dispatchGasConcentrationsToDEM())
    {
        fatalErrorInFunction << "dispatchGasConcentrationsToDEM() failed.\n";
        Plus::processor::abort(0);
    }

    sendDataTimer().end();
    return true;
}

Foam::tmp<Foam::volScalarField>
multiReactiveThermalMomentumSphereUnresolvedCouplingSystem::massSource
(
    const word& specieName
) const
{
    const auto& names = reactionInteraction_.gasSpeciesNames();

    forAll(names, g)
    {
        if (names[g] == specieName)
        {
            return Foam::tmp<Foam::volScalarField>
            (
                reactionInteraction_.speciesSu()[g]
            );
        }
    }

    FatalErrorInFunction
        << "massSource() requested for unknown species '" << specieName
        << "'. Registered reactive gas species are: " << names
        << Foam::exit(Foam::FatalError);

    return Foam::tmp<Foam::volScalarField>(nullptr);
}

} // pFlow::coupling
