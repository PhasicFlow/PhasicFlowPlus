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

template<typename DragClosureType>
pFlow::coupling::grainDrag<DragClosureType>::grainDrag
(
    const unresolvedCouplingSystem& uCS, 
    const porosity& 				prsty
)
:
    drag(uCS, prsty),
    dragClosure_(this->dict()),
    courseGrainFactor_
    (
        static_cast<const momentumGrainUnresolvedCouplingSystem&>(uCS).courseGrainFactor()
    )
{}

template<typename DragClosureType>
void pFlow::coupling::grainDrag<DragClosureType>::calculateDragForce
(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    diameter,
    const distributionBase&         cellDistribution,    
    Plus::realx3ProcCMField&        particleForce
)
{
    setSuSpToZero();

    const auto& parCellInd =  this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    size_t numPar = parCellInd.size();
    const auto& alpha = this->alpha();
    
    auto fluidVel = fluidVelocity.fieldSpan();
    auto solidVel = parVelocity.fieldSpan();

    auto& Su = this->Su();
    auto& Sp = this->Sp();

    // gets pressure gradient 
    auto pGradPtr = this->pressureGradient(rho);
    const auto& pGrad = pGradPtr();

    #pragma omp parallel for schedule (dynamic)
    for(size_t parIndx=0; parIndx<numPar; parIndx++)
    {
        auto cellIndx = parCellInd[parIndx];

        if(cellIndx < 0 ) continue;

        auto rhoi = rho[cellIndx];
        auto mui = nu[cellIndx]* rhoi;
        auto ef = alpha[cellIndx];
        auto dp = diameter[parIndx];
        auto cgf = courseGrainFactor_[parIndx];
        auto dps = dp/cgf;
        auto vp =  Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);

        Foam::vector up{solidVel[parIndx].x(), solidVel[parIndx].y(),solidVel[parIndx].z()};

        Foam::vector ur = fluidVel[parIndx]-up;

        Foam::scalar Res = ef * rhoi * Foam::mag(ur) * dps /mui ;

        Foam::scalar sp = 3 * Foam::pow(cgf,3) * Foam::constant::mathematical::pi * 
                    mui * ef * dps * dragClosure_.dimlessDrag(Res, ef);
            
        Foam::vector pf = static_cast<real>(sp)*ur - vp*pGrad[cellIndx];

        particleForce[parIndx] += realx3(pf.x(), pf.y(), pf.z());
        
        
        cellDistribution.distributeValue_OMP(parIndx, cellIndx, Su, -(sp*up));
        cellDistribution.distributeValue_OMP(parIndx, cellIndx, Sp,   sp);
        
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
