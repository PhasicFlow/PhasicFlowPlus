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
#include "celloPolydisperseSphereDrag.hpp"
#include "distributionBase.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"


pFlow::coupling::celloPolydisperseSphereDrag::celloPolydisperseSphereDrag
(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty
)
:
    polydisperseDragBase(uCS, prsty),
    dragClosure_(this->dict())
{}

void pFlow::coupling::celloPolydisperseSphereDrag::calculateDragForce
(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,
    Plus::realx3ProcCMField&        particleForce
)
{
    setSuSpToZero();

    const auto& parCellInd = this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    const auto& alpha = this->alpha();

    auto fluidVel = fluidVelocity.fieldSpan();
    auto solidVel = parVelocity.fieldSpan();

    auto& Su = this->Su();
    auto& Sp = this->Sp();

    // gets pressure gradient
    auto pGradPtr = this->pressureGradient(rho);
    const auto& pGrad = pGradPtr();

    // sauter Diameter call
    calculateSauterDiameter(diameter, cellDistribution);
    const auto& dBar = this->averageDiameter();
    const auto& sumD3 = this->sumD3();
    const auto& sumD4 = this->sumD4();
    

    const size_t numPar = parCellInd.size();

   

    #pragma omp parallel for schedule(dynamic)
    for(size_t parIndx=0; parIndx<numPar; parIndx++)
    {
        const auto cellIndx = parCellInd[parIndx];

        if(cellIndx < 0 ) continue;

        const Foam::scalar rhoi = rho[cellIndx];
        const Foam::scalar mui = nu[cellIndx]*rhoi;
        const Foam::scalar ef = Foam::max(alpha[cellIndx], Foam::SMALL);
        const Foam::scalar dp = diameter[parIndx];
        const Foam::scalar dRef = Foam::max(dBar[cellIndx], Foam::SMALL);

        const Foam::scalar vp =  Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);

        Foam::vector up{solidVel[parIndx].x(), solidVel[parIndx].y(), solidVel[parIndx].z()};
        Foam::vector ur = fluidVel[parIndx]-up;

        // Step 2: evaluate Cello closure with reference diameter dBar
        const Foam::scalar ReBar = ef * rhoi * Foam::mag(ur) * dRef /mui;
        const Foam::scalar fBar = dragClosure_.dimlessDrag(ReBar, ef);

        // Step 3: beta_i
        const Foam::scalar yi = dp/dRef;
        const Foam::scalar sumXkYk =
            sumD3[cellIndx] > Foam::SMALL
          ? (sumD4[cellIndx]/(sumD3[cellIndx]*dRef))
          : Foam::scalar(1);

        const Foam::scalar beta = yi +
            ((1.0-ef)/ef)*((1.0-ef-0.27)/(1.0-0.27))*
            ((Foam::sqr(yi)-yi)/Foam::max(sumXkYk, Foam::SMALL));


        const Foam::scalar sp = beta * 3.0 * Foam::constant::mathematical::pi * mui * ef * dp * fBar;

        const Foam::vector pf = static_cast<real>(sp)*ur - vp*pGrad[cellIndx];

        particleForce[parIndx] += realx3(pf.x(), pf.y(), pf.z());

        cellDistribution.distributeValue_OMP(parIndx, cellIndx, Su, -(sp*up));
        cellDistribution.distributeValue_OMP(parIndx, cellIndx, Sp, sp);
    }

    const auto& Vcells = this->mesh().V();

    forAll(Vcells, i)
    {
        Su[i] /= Vcells[i];
        Sp[i] /= Vcells[i];
    }

    cellDistribution.smoothenField(Sp);
    cellDistribution.smoothenField(Su);

    Sp.correctBoundaryConditions();
    Su.correctBoundaryConditions();
}

