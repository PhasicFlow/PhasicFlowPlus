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

#include "thermalMomentumSphereUnresolvedCouplingSystem.hpp"
#include <algorithm>   

namespace pFlow
{
namespace coupling
{

//----------------------------- private methods -------------------------------

bool thermalMomentumSphereUnresolvedCouplingSystem::collectFluidHeatSource()
{
    auto& comm = parMapping().realScatteredComm();

    auto allHeatConv = pDEMSystem().parFluidHeatSourceConv();

    if (Plus::processor::isMaster() && allHeatConv.size() > 0)
    {
         std::fill(
             allHeatConv.data(), 
             allHeatConv.data() + allHeatConv.size(), 
             real(0));
    }

    auto thisHeatConv = makeSpan(fluidHeatSourceConv_);
    if (!comm.collectSum(thisHeatConv, allHeatConv)) 
    {
        return false;
    }

    auto allHeatRad = pDEMSystem().parFluidHeatSourceRad();

    if (Plus::processor::isMaster() && allHeatRad.size() > 0)
    {
         std::fill(
             allHeatRad.data(), 
             allHeatRad.data() + allHeatRad.size(), 
             real(0));
    }

    auto thisHeatRad = makeSpan(fluidHeatSourceRad_);
    if (!comm.collectSum(thisHeatRad, allHeatRad)) 
    {
        return false;
    }

    return true;
}

bool thermalMomentumSphereUnresolvedCouplingSystem::collectFluidProperties()
{
    auto& comm = parMapping().realScatteredComm();

    auto allFluidKappa = pDEMSystem().parFluidKappa();

    if (Plus::processor::isMaster() && allFluidKappa.size() > 0)
    {
         std::fill(
             allFluidKappa.data(), 
             allFluidKappa.data() + allFluidKappa.size(), 
             real(0));
    }

    auto thisFluidKappa = makeSpan(fluidKappa_);
    if (!comm.collectSum(thisFluidKappa, allFluidKappa)) 
    {
        return false;
    }

    auto allFluidAlpha = pDEMSystem().parFluidAlpha();

    if (Plus::processor::isMaster() && allFluidAlpha.size() > 0)
    {
         std::fill(
             allFluidAlpha.data(), 
             allFluidAlpha.data() + allFluidAlpha.size(), 
             real(0));
    }

    auto thisFluidAlpha = makeSpan(fluidAlpha_);
    if (!comm.collectSum(thisFluidAlpha, allFluidAlpha)) 
    {
        return false;
    }

    if (Plus::processor::isMaster() && allFluidKappa.size() > 0)
    {
        const size_t N  = allFluidKappa.size();
        size_t badKappa = 0;
        size_t badAlpha = 0;

        for (size_t i = 0; i < N; ++i)
        {
            if (allFluidKappa[i] <= real(0))
            {
                ++badKappa;
                allFluidKappa[i] = real(1e-8);  
            }

            const real a = allFluidAlpha[i];
            if (a < real(0) || a > real(1))
            {
                ++badAlpha;
                allFluidAlpha[i] = std::max(
                    std::min(a, real(1)), 
                    real(0));
            }
        }

        // Per-step diagnostic: flags a single step where clamping is
        // already widespread (>0.1% of particles).
        const real fracKappa = 
            static_cast<real>(badKappa) / static_cast<real>(N);
            
        const real fracAlpha = 
            static_cast<real>(badAlpha) / static_cast<real>(N);

        if (fracKappa > real(0.001))
        {
            pFlow::output
                << "[WARNING thermalCoupling] " << badKappa
                << " particles (" << 100.0 * fracKappa
                << "%) have kappa_fluid <= 0 after collection. "
                << "Clamped to 1e-8 W/(m.K).\n";
        }

        if (fracAlpha > real(0.001))
        {
            pFlow::output
                << "[WARNING thermalCoupling] " << badAlpha
                << " particles (" << 100.0 * fracAlpha
                << "%) have alpha outside [0,1] after collection. Clamped.\n";
        }

        // Cumulative diagnostic: catches a persistent low-level issue
        // that never crosses the per-step threshold above.
        cumulativeBadKappaCount_ += badKappa;
        cumulativeBadAlphaCount_ += badAlpha;

        const uint64 cumulativeTotal = 
            cumulativeBadKappaCount_ + cumulativeBadAlphaCount_;

        if (cumulativeTotal >= nextKappaAlphaReportMilestone_)
        {
            pFlow::output
                << "[WARNING thermalCoupling] Cumulative clamping since run "
                << "start: " << cumulativeBadKappaCount_ 
                << " particle-steps with kappa_fluid <= 0, "
                << cumulativeBadAlphaCount_ 
                << " particle-steps with alpha outside [0,1].\n";

            // Next reminder at double the current milestone, so reminders
            // become progressively less frequent for a long-running case.
            nextKappaAlphaReportMilestone_ = cumulativeTotal * 2;
        }
    }

    return true;
}

void thermalMomentumSphereUnresolvedCouplingSystem::sendFluidHeatSourceToDEM()
{
    collectFluidHeatSource();

    if (!pDEMSystem().sendFluidHeatSourcesToDEM())
    {
        fatalErrorInFunction << "sendFluidHeatSourcesToDEM failed.\n";
        Plus::processor::abort(0);
    }
}

void thermalMomentumSphereUnresolvedCouplingSystem::sendFluidPropertiesToDEM()
{
    const auto& kappa   = this->cMesh().mesh()
                              .lookupObject<Foam::volScalarField>("kappa");
    const auto& alpha   = this->alpha();
    const auto& cellIDs = this->parCellIndex();
    const Foam::label nCells = static_cast<Foam::label>(kappa.size());

    size_t invalidCellCount = 0;

    for (size_t i = 0; i < fluidKappa_.size(); ++i)
    {
        const Foam::label cellI = cellIDs[i];

        if (cellI >= 0 && cellI < nCells)
        {
            fluidKappa_[i] = static_cast<real>(kappa [cellI]);
            fluidAlpha_[i] = static_cast<real>(alpha [cellI]);
        }
        else
        {
            fluidKappa_[i] = real(0);
            fluidAlpha_[i] = real(0);
            ++invalidCellCount;
        }
    }

    if (fluidKappa_.size() > 0)
    {
        // Per-step diagnostic.
        const real invalidFrac =
            static_cast<real>(invalidCellCount) / 
            static_cast<real>(fluidKappa_.size());

        if (invalidFrac > real(0.05))
        {
            pFlow::output
                << "[WARNING thermalCoupling] " << invalidCellCount
                << " of " << fluidKappa_.size()
                << " particles (" << 100.0 * invalidFrac
                << "%) have cellI < 0 in sendFluidPropertiesToDEM().\n";
        }

        // Cumulative diagnostic: same pattern as collectFluidProperties()
        // above -- catches a persistent trickle of invalid-cell particles
        // (fluidKappa zeroed, excluded from PFP averaging) below the
        // per-step threshold.
        cumulativeInvalidCellCount_ += invalidCellCount;

        if (cumulativeInvalidCellCount_ >= nextInvalidCellReportMilestone_)
        {
            pFlow::output
                << "[WARNING thermalCoupling] Cumulative invalid-cell count "
                << "since run start: " << cumulativeInvalidCellCount_
                << " particle-steps with cellI < 0 in "
                << "sendFluidPropertiesToDEM() (fluidKappa/fluidAlpha zeroed "
                << "and excluded from PFP averaging for each).\n";

            nextInvalidCellReportMilestone_ = cumulativeInvalidCellCount_ * 2;
        }
    }

    collectFluidProperties();

    if (!pDEMSystem().sendFluidPropertiesToDEM())
    {
        fatalErrorInFunction << "sendFluidPropertiesToDEM() failed.\n";
        Plus::processor::abort(0);
    }
}

//----------------------------- protected methods -----------------------------

bool thermalMomentumSphereUnresolvedCouplingSystem::distributeParticleFields()
{
    if (!momentumSphereUnresolvedCouplingSystem::distributeParticleFields())
    {
        return false;
    }

    auto& comm = parMapping().realScatteredComm();

    auto allTemp  = pDEMSystem().temperature();
    auto thisTemp = makeSpan(particleTemperature_);
    if (!comm.distribute(allTemp, thisTemp)) 
    {
        return false;
    }

    auto allEmis  = pDEMSystem().emissivity();
    auto thisEmis = makeSpan(emissivity_);
    if (!comm.distribute(allEmis, thisEmis)) 
    {
        return false;
    }

    // Radiation's neighbourhood-sum fields are only distributed when
    // radiation is actually enabled. The DEM side keeps these
    // correctly zero-filled (not empty) even when disabled -- see
    // thermalSphereDEMSystem::ensureRadiationHostMemory() -- so
    // skipping this block is purely an efficiency choice (no MPI
    // traffic for a mechanism nobody is using this run), not a
    // correctness requirement on its own; local buffers are kept at
    // their well-defined zero in the disabled case below regardless.
    if (pDEMSystem().hasRadiation())
    {
        auto allRadSum  = pDEMSystem().radSumTemp();
        auto thisRadSum = makeSpan(radSumTemp_);
        if (!comm.distribute(allRadSum, thisRadSum)) 
        {
            return false;
        }

        auto allNum = pDEMSystem().radNumPrt();

        std::vector<real> dNumPrtR;
        if (Plus::processor::isMaster())
        {
            dNumPrtR.resize(allNum.size());
            for (size_t i = 0; i < allNum.size(); ++i)
            {
                dNumPrtR[i] = static_cast<real>(allNum[i]);
            }
        }

        std::vector<real> pNumPrtR(radNumPrt_.size(), real(0));
        pFlow::span<real> dCountR(dNumPrtR.data(), dNumPrtR.size());
        pFlow::span<real> pCountR(pNumPrtR.data(), pNumPrtR.size());
        
        if (!comm.distribute(dCountR, pCountR)) 
        {
            return false;
        }

        for (size_t i = 0; i < pNumPrtR.size(); ++i)
        {
            radNumPrt_[i] = static_cast<pFlow::uint32>(pNumPrtR[i]);
        }
    }
    else
    {
        std::fill(radSumTemp_.begin(), radSumTemp_.end(), real(0));
        std::fill(radNumPrt_.begin(), radNumPrt_.end(), pFlow::uint32(0));
    }

    return true;
}

//----------------------------- constructors ----------------------------------

thermalMomentumSphereUnresolvedCouplingSystem::
thermalMomentumSphereUnresolvedCouplingSystem(
    word            shapeTypeName,
    word            couplingSystemType,
    Foam::fvMesh&   mesh,
    int             argc,
    char*           argv[])
:
    momentumSphereUnresolvedCouplingSystem(
        shapeTypeName, 
        couplingSystemType, 
        mesh, 
        argc, 
        argv),
    heatInteraction_(
        *this, 
        porosityCoupling()),
    particleTemperature_(
        "temperature",    
        this->parMapping().centerMass()),
    fluidHeatSourceConv_(
        "heatSourceConv", 
        this->parMapping().centerMass()),
    fluidHeatSourceRad_(
        "heatSourceRad",  
        this->parMapping().centerMass()),
    emissivity_(
        "emissivity",     
        this->parMapping().centerMass()),
    radSumTemp_(
        "radSumTemp",     
        this->parMapping().centerMass()),
    radNumPrt_(
        "radNumPrt",      
        this->parMapping().centerMass()),
    fluidKappa_(
        "fluidKappa",     
        this->parMapping().centerMass()),
    fluidAlpha_(
        "fluidAlpha",     
        this->parMapping().centerMass())
{
    requiresDistribution_ = heatInteraction_.requireCellDistribution();
}

//---------------------------- public methods ---------------------------------

void thermalMomentumSphereUnresolvedCouplingSystem::calculatePorosity()
{
    momentumSphereUnresolvedCouplingSystem::calculatePorosity();

    // The base class call above only refreshes distribution weights when
    // ITS OWN requirement (porosity + momentum) calls for it. If heat is
    // the only side that needs them, refresh here too, before
    // calculateHeatCoupling() runs later this step.
    if (heatInteraction_.requireCellDistribution() &&
        !momentumSphereUnresolvedCouplingSystem::requireCellDistribution())
    {
        this->updateDistributionWeights();
    }
}

// calculateMomentumCoupling() is no longer overridden here - it is
// inherited unchanged from momentumSphereUnresolvedCouplingSystem, which
// now owns the same momentumInteraction_ this class used to duplicate.

void thermalMomentumSphereUnresolvedCouplingSystem::calculateHeatCoupling()
{
    const auto& U  = this->cMesh().mesh()
                          .lookupObject<Foam::volVectorField>("U");
    const auto& vp = this->particleVelocity();

    // Momentum's already-computed averaging, used as a fallback when heat
    // doesn't own its own (similarToMomentum) - see heatInteraction_'s own
    // calculateCoupling() for which one actually gets used.
    const auto& fluidVel = momentumInteractionCoupling().fluidVelAveraging();
    const auto& parVel   = momentumInteractionCoupling().solidVelAveraging();

    const auto& dp = this->particleDiameter();
    const auto& Tp = this->particleTemperature();

    // Zero local heat source buffers before new calculation
    std::fill(fluidHeatSourceConv_.begin(), fluidHeatSourceConv_.end(), 0.0);
    std::fill(fluidHeatSourceRad_ .begin(), fluidHeatSourceRad_ .end(), 0.0);

    heatInteraction_.calculateCoupling(
        U,
        vp,
        fluidVel,
        parVel,
        dp,
        Tp,
        fluidHeatSourceConv_,   
        emissivity_,
        radSumTemp_,
        radNumPrt_,
        fluidHeatSourceRad_);
}

// calculateMassCoupling(), Sp() and Su() are no longer overridden here
// either, for the same reason as calculateMomentumCoupling() above -
// they are inherited unchanged from momentumSphereUnresolvedCouplingSystem.

Foam::tmp<Foam::volScalarField>
thermalMomentumSphereUnresolvedCouplingSystem::heatSource() const
{
    return Foam::tmp<Foam::volScalarField>(heatInteraction_.Sh());
}

Foam::tmp<Foam::volVectorField>
thermalMomentumSphereUnresolvedCouplingSystem::Us() const
{
    // Return the Eulerian solid-phase velocity field if registered
    if (this->cMesh().mesh().foundObject<Foam::volVectorField>("Us"))
    {
        return Foam::tmp<Foam::volVectorField>(
            this->cMesh().mesh().lookupObject<Foam::volVectorField>("Us"));
    }

    // Graceful fallback: zero velocity
    return momentumSphereUnresolvedCouplingSystem::Us();
}

bool thermalMomentumSphereUnresolvedCouplingSystem::sendDataToDEM(
    real t, 
    real dt)
{
    if (!momentumSphereUnresolvedCouplingSystem::sendDataToDEM(t, dt))
    {
        return false;
    }

    sendDataTimer().start();

    sendFluidHeatSourceToDEM();
    sendFluidPropertiesToDEM();

    sendDataTimer().end();

    return true;
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling
} // pFlow
