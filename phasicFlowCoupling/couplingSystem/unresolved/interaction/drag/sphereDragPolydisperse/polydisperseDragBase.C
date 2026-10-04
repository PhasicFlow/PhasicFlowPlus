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
 
#include "polydisperseDragBase.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "distributionBase.hpp"
 
pFlow::coupling::polydisperseDragBase::polydisperseDragBase
(
     const unresolvedCouplingSystem& uCS, 
     const porosity& 				prsty
)
 :
     drag(uCS, prsty),
     diameterClass_("diameterClass", uCS.centerMass()),
     writeAverageDiameter_
     (
         this->dict().template lookupOrDefault<Foam::Switch>
         (
             "writeAverageDiameter",
             Foam::Switch(false)
         )
     ),
     averageDiameter_
     (
         Foam::IOobject
         (
             "averageDiameter",
             mesh().time().timeName(),
             mesh(),
             Foam::IOobject::NO_READ,
             (writeAverageDiameter_?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
         ),
         mesh(),
         Foam::dimensionedScalar(Foam::dimLength, Foam::scalar(0))
     ),
    sumD2_
    (
    	Foam::IOobject
    	(
        "sumD2",
        mesh().time().timeName(),
        mesh(),
        Foam::IOobject::NO_READ,
        Foam::IOobject::NO_WRITE,
        false
    	),
    	mesh(),
    	Foam::dimensionedScalar(Foam::pow(Foam::dimLength,2), Foam::scalar(0))
     ),
    sumD3_
    (
    	Foam::IOobject
    	(
        "sumD3",
        mesh().time().timeName(),
        mesh(),
        Foam::IOobject::NO_READ,
        Foam::IOobject::NO_WRITE,
        false
    	),
    	mesh(),
    	Foam::dimensionedScalar(Foam::pow(Foam::dimLength,3), Foam::scalar(0))
     ),
    sumD4_
     (
    	Foam::IOobject
    	(
        "sumD4",
        mesh().time().timeName(),
        mesh(),
        Foam::IOobject::NO_READ,
        Foam::IOobject::NO_WRITE,
        false
    	),
    	mesh(),
    	Foam::dimensionedScalar(Foam::pow(Foam::dimLength,4), Foam::scalar(0))
      )

 {}


void pFlow::coupling::polydisperseDragBase::calculateSauterDiameter
(
    const Plus::realProcCMField& diameter,
    const distributionBase&      cellDistribution
)
{
    auto& d32 = averageDiameter_;
    auto& d2Sum = sumD2_;
    auto& d3Sum = sumD3_;
    auto& d4Sum = sumD4_;
    
    d32 = Foam::dimensionedScalar(Foam::dimLength, Foam::scalar(0));
    d2Sum = dimensionedScalar(pow(dimLength,2), scalar(0));
    d3Sum = dimensionedScalar(pow(dimLength,3), scalar(0));
    d4Sum = dimensionedScalar(pow(dimLength,4), scalar(0));
    
    

    const auto& parCellInd = this->parCellIndex();
    const auto numPar = parCellInd.size();

    #pragma omp parallel for schedule(dynamic)
    for(size_t parIndx = 0; parIndx < numPar; parIndx++)
    {
        const auto cellIndx = parCellInd[parIndx];

        if(cellIndx < 0) continue;

        const Foam::scalar dp = diameter[parIndx];
        const Foam::scalar dp2 = Foam::sqr(dp);
        const Foam::scalar dp3 = dp2*dp;
        const Foam::scalar dp4 = dp2*dp2;
	
	cellDistribution.distributeValue_OMP(parIndx, cellIndx, d2Sum, dp2);
	cellDistribution.distributeValue_OMP(parIndx, cellIndx, d4Sum, dp4);
        cellDistribution.distributeValue_OMP(parIndx, cellIndx, d3Sum, dp3);
        
        
    }

    cellDistribution.smoothenField(d2Sum);
    cellDistribution.smoothenField(d4Sum);
    cellDistribution.smoothenField(d3Sum);
    

    

    forAll(d32, cellI)
    {
        if(d2Sum[cellI] > Foam::SMALL)
        {
            d32[cellI] = d3Sum[cellI]/d2Sum[cellI];
        }
        else
        {
            d32[cellI] = Foam::scalar(0);
        }
    }

    d32.correctBoundaryConditions();

}
