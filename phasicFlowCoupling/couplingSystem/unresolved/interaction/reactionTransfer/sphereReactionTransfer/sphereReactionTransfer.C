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

#include "sphereReactionTransfer.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

//----------------------------- constructors ----------------------------------

sphereReactionTransfer::sphereReactionTransfer
(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty
)
:
    reactionTransfer(uCS, prsty)
{
    // Field allocation is handled by the reactionTransfer base constructor
}


//---------------------------- public methods ---------------------------------

void sphereReactionTransfer::calculateReactionTransfer
(
    const std::vector<real>&     gasMassSource,
    const std::vector<real>&     gasMassSourceSp,
    const distributionBase&      cellDistribution
)
{
    // Zero all source fields before the accumulation
    setSuSpToZero();

    const auto& parCellInd = this->parCellIndex();
    const auto& Vcells     = this->mesh().V();

    const size_t numPar = parCellInd.size();
    const size_t N_gas  = this->gasSpeciesNames().size();

    auto& suSpec = this->speciesSu();
    auto& spSpec = this->speciesSp();

    // Distribute the particle sources to the CFD cells. Thread safety relies
    // on the OpenMP atomics inside distributeValue_OMP.
    #pragma omp parallel for schedule(dynamic)
    for (size_t parIndx = 0; parIndx < numPar; ++parIndx)
    {
        const auto cellIndx = parCellInd[parIndx];

        // Skip particles outside this processor's local fluid mesh domain
        if (cellIndx < 0) continue;

        // Distribute gas species mass sources per species [kg/s]
        for (size_t g = 0; g < N_gas; ++g)
        {
            const Foam::scalar m_dot_g = gasMassSource[parIndx * N_gas + g];

            cellDistribution.distributeValue_OMP
            (
                parIndx, cellIndx, suSpec[g], m_dot_g
            );

            const Foam::scalar sp_g = gasMassSourceSp[parIndx * N_gas + g];

            cellDistribution.distributeValue_OMP
            (
                parIndx, cellIndx, spSpec[g], sp_g
            );
        }
    }

    // Convert the integrated totals to volumetric source terms
    forAll(Vcells, i)
    {
        for (size_t g = 0; g < N_gas; ++g)
        {
            suSpec[g][i] /= Vcells[i];   // [kg/s]           -> [kg/(m3.s)]
            spSpec[g][i] /= Vcells[i];   // [kg.m3/(kmol.s)] -> [kg/(kmol.s)]
        }
    }

    // Smooth the fields and synchronise them across the MPI boundaries
    for (size_t g = 0; g < N_gas; ++g)
    {
        cellDistribution.smoothenField(suSpec[g]);
        suSpec[g].correctBoundaryConditions();

        cellDistribution.smoothenField(spSpec[g]);
        spSpec[g].correctBoundaryConditions();
    }
}

} // namespace pFlow::coupling
