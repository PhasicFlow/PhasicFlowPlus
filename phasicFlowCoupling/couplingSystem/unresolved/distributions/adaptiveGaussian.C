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

#define distributed1

#include "adaptiveGaussian.hpp"

pFlow::coupling::adaptiveGaussian::adaptiveGaussian
(
    Foam::dictionary dict, 
    const couplingMesh &cMesh, 
    const Plus::centerMassField &centerMass
)
:
    distribution(dict, cMesh, centerMass),
    distLength_(lookupOrDefaultDict(dict, "distLength", static_cast<Foam::scalar>(1.0))),
    maxLayers_(lookupOrDefaultDict(dict, "maxLayers", static_cast<Foam::label>(2)))
{
    constructLists(distLength_, maxLayers_);
}

void pFlow::coupling::adaptiveGaussian::updateWeights
(
	const Plus::procCMField<Foam::label> & parCellIndex,
	const Plus::procCMField<real> & parDiameter
)
{
	const auto& centerMass = weights_.centerMass();
	const size_t numPar = centerMass.size();
	const Foam::scalarField& cellV = mesh_.cellVolumes();
	const Foam::vectorField& cellC = mesh_.cellCentres(); 

  //#pragma omp parallel for
	for(size_t i=0; i<numPar; i++)
	{
		const Foam::label targetCellId = parCellIndex[i];
		auto& parWeights = weights_[i];
		
		parWeights.clear();
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const auto& neighbors = neighborList_[targetCellId];

		const Foam::scalar dp = parDiameter[i];
		const Foam::scalar dcell = Foam::pow(cellV[targetCellId], 0.333333);

		const realx3& cp_i = centerMass[i];
    const Foam::vector cp{cp_i.x(), cp_i.y(), cp_i.z()};
    
    const Foam::scalar dx_dp = dcell/dp;
    // work like PIC
    if(dx_dp>10.0)
    {
        parWeights.push_back({targetCellId,1.0});
        continue;      
    }
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
      if( f>1.0e-4)
      {
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
    if( Foam::scalar f = Foam::exp(-0.5*dist2/std2); f>1.0e-5)
    {
      parWeights.push_back({cellId,f});
      pSubTotal += f;
    }
  }
#endif
		pSubTotal = Foam::max(pSubTotal, static_cast<Foam::scalar>(1.0e-10));
		for(auto& [cellid, w]:parWeights) w /= pSubTotal;     
		
	}
	
}
