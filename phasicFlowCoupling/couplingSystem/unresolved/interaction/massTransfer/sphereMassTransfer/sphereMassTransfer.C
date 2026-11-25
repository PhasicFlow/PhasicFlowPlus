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

#include "sphereMassTransfer.hpp"

namespace pFlow
{

namespace coupling
{

//----------------------------- constructors ----------------------------------

template<typename SherwoodClosureType>
sphereMassTransfer<SherwoodClosureType>::sphereMassTransfer(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty)
:
    massTransfer(uCS, prsty),
    sherwoodClosure_(this->dict())
{}

//---------------------------- public methods ---------------------------------

template<typename SherwoodClosureType>
void sphereMassTransfer<SherwoodClosureType>::calculateMassTransferCoeff(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    diameter,
    std::vector<real>&              massTransferCoeff)
{
    const auto& parCellInd = this->parCellIndex();
    const size_t numPar    = parCellInd.size();

    massTransferCoeff.assign(numPar, real(0));

    const auto& nu =
        this->mesh().template lookupObject<Foam::volScalarField>("nu");

    const auto& tp = this->mesh().template lookupObject<Foam::IOdictionary>(
        "transportProperties");

    const Foam::scalar Sc =
        Foam::dimensionedScalar("Sc", Foam::dimless, tp).value();

    if (Sc <= 0)
    {
        FatalErrorInFunction
            << "Sc in transportProperties must be positive, not " << Sc
            << Foam::exit(Foam::FatalError);
    }

    const auto& alpha = this->alpha();

    auto fluidVel = fluidVelocity.fieldSpan();
    auto solidVel = parVelocity.fieldSpan();

    // Lower bounds, read outside the loop
    const Foam::scalar nuMin =
        this->dict().lookupOrDefault<Foam::scalar>("nuMin", 1e-15);

    const Foam::scalar dpMin =
        this->dict().lookupOrDefault<Foam::scalar>("dpMin", 1e-15);

    const Foam::scalar alphaMin = this->Porosity().alphaMin();

    #pragma omp parallel for schedule(dynamic)
    for (size_t parIndx = 0; parIndx < numPar; ++parIndx)
    {
        const auto cellIndx = parCellInd[parIndx];

        // Particles outside this processor's local fluid mesh domain
        if (cellIndx < 0) continue;

        const Foam::scalar nui = Foam::max(nu[cellIndx], nuMin);
        const Foam::scalar dp  = Foam::max(diameter[parIndx], dpMin);

        const Foam::scalar ef =
            Foam::min(Foam::max(alpha[cellIndx], alphaMin), 1.0);

        // Relative slip velocity
        const Foam::vector up
        {
            solidVel[parIndx].x(),
            solidVel[parIndx].y(),
            solidVel[parIndx].z()
        };
        const Foam::vector ur = fluidVel[parIndx] - up;

        const Foam::scalar Re = ef * Foam::mag(ur) * dp / nui;
        const Foam::scalar Sh = sherwoodClosure_(Re, Sc, ef);

        // k_c = Sh*D/d_p with D = nu/Sc
        massTransferCoeff[parIndx] = static_cast<real>(Sh * nui / (Sc * dp));
    }
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling

} // pFlow
