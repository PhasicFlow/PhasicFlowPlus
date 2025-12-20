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

#include "adaptiveGaussian.hpp"
#include "couplingMesh.hpp"
#include "streams.hpp"


pFlow::coupling::adaptiveGaussian::adaptiveGaussian
(
    const Foam::dictionary& 	 parrentDict, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    distribution(parrentDict, cMesh, centerMass)
{
    
    auto dict = parrentDict.subOrEmptyDict("adaptiveGaussianInfo");

    maxLayers_ = dict.getOrDefault(
            "maxLayers", 
            static_cast<Foam::label>(1));

    smoothingFactor_ =  dict.getOrDefault
    (
        "smoothingFactor", 
        static_cast<Foam::scalar>(1.0)
    );

    REPORT(1) << "Smoothing factor for adaptiveGaussian distribution is: "
        <<  Yellow_Text(smoothingFactor_) << END_REPORT;
    
    
}

void pFlow::coupling::adaptiveGaussian::checkForListsConstructed()
{
    if( !listsConstructed_ )
    {
        REPORT(0)<<"Constructing neighbor lists for adaptiveGaussian distribution..."<<END_REPORT;

        constructLists(1.0, maxLayers_);
        listsConstructed_ = true;
    }
}

void pFlow::coupling::adaptiveGaussian::updateWeights
(
    const Plus::procCMField<real> & parDiameter
)
{
    
    checkForListsConstructed();

    const auto& parCellIndex = cMesh().parCellIndex();
    auto& weights = this->weights();
    const auto& centerMass = weights.centerMass();
    const size_t numPar = centerMass.size();
    const Foam::scalarField& cellV = mesh().cellVolumes();
    const Foam::vectorField& cellC = mesh().cellCentres(); 

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        const Foam::label targetCellId = parCellIndex[i];
        auto& parWeights = weights[i];
        
        parWeights.clear();
        if( targetCellId < 0 )continue;
        
        // get all the neighbors of cell 
        const auto& neighbors = neighborList_[targetCellId];

        const Foam::scalar dp = parDiameter[i];
        const Foam::scalar dcell = Foam::pow(cellV[targetCellId], 0.333333);

        const realx3& cp_i = centerMass[i];
        const Foam::vector cp{cp_i.x(), cp_i.y(), cp_i.z()};
        
        const Foam::scalar dx_dp = dcell/dp;
        
        // work like PCM
        if(dx_dp>7.0)
        {
            parWeights.push_back({targetCellId,1.0});
            continue;      
        }

        const auto std2 = Foam::pow(stdDeviation(dx_dp,dcell),2);

        Foam::scalar pSubTotal = 0;
        for(auto cellId:neighbors)
        {
            Foam::vector xx = cellC[cellId] - cp;
            Foam::scalar dist2 = Foam::dot(xx,xx);
            if( Foam::scalar f = Foam::exp(-0.5*dist2/std2); f>1.0e-3)
            {
                parWeights.push_back({cellId,f});
                pSubTotal += f;
            }
        }

        pSubTotal = Foam::max(pSubTotal, static_cast<Foam::scalar>(1.0e-10));
        for(auto& [cellid, w]:parWeights) w /= pSubTotal;     
        
    }
    
}


/* this is for future references 
#ifdef distributed1
    
        Foam::scalar pSubTotal = 0;
        for(auto cellId:neighbors)
        {
            auto dcell_D = Foam::pow(cellV[cellId], 0.333333);
            auto dx_dp_D = dcell_D/ dp ;
            auto std2_D = Foam::pow(stdDeviation(dx_dp_D,dcell_D),2);
            
            Foam::vector xx = cellC[cellId] - cp;
            Foam::scalar dist2 = Foam::dot(xx,xx);
            Foam::scalar f= Foam::exp(-0.5*dist2/std2_D);
            if( f>1.0e-3)
            {
                //selected ++;
                parWeights.push_back({cellId,f});
                pSubTotal += f;
            }
        }
    
#else
    const auto std2 = Foam::pow(stdDeviation(dx_dp,dcell),2);

    Foam::scalar pSubTotal = 0;
    for(auto cellId:neighbors)
    {
        Foam::vector xx = cellC[cellId] - cp;
        Foam::scalar dist2 = Foam::dot(xx,xx);
        if( Foam::scalar f = Foam::exp(-0.5*dist2/std2); f>1.0e-3)
        {
            parWeights.push_back({cellId,f});
            pSubTotal += f;
        }
    }
#endif
*/