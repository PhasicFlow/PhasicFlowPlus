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

// from OpenFOAM
#include "GaussianIntegral.hpp"
#include "couplingMesh.hpp"
#include "schedule.hpp"


pFlow::coupling::GaussianIntegral::GaussianIntegral
(
	Foam::dictionary 		dict, 
	const couplingMesh& 	cMesh,
	const Plus::centerMassField& centerMass
)
:
	distribution(dict, cMesh, centerMass),
	maxLayers_(lookupOrDefaultDict(dict, "maxLayers", static_cast<Foam::label>(1)))
{
	constructLists(1.0, maxLayers_);
}


void pFlow::coupling::GaussianIntegral::updateWeights
(
	const Plus::procCMField<Foam::label> & parCellIndex,
	const Plus::procCMField<real> & parDiameter
)
{
	const auto& centerMass = weights_.centerMass();
	const size_t numPar = centerMass.size();
	const Foam::scalarField& cellV = mesh_.cellVolumes();
	const Foam::vectorField& cellC = mesh_.cellCentres(); 

	#pragma omp parallel for schedule (dynamic)
	for(size_t i=0; i<numPar; i++)
	{
		const Foam::label targetCellId = parCellIndex[i];
		auto& parWeights = weights_[i];
		
		parWeights.clear();
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const auto& neighbors = neighborList_[targetCellId];

		const Foam::scalar rp = parDiameter[i]/2;
		const Foam::scalar vp = 4*Pi/3 * Foam::pow(rp,3);
		const realx3& cp_i = centerMass[i];
		const Foam::vector cp{cp_i.x(), cp_i.y(), cp_i.z()};

		Foam::scalar pSubTotal = 0;
		for(auto cellId:neighbors)
		{

			const auto vc = cellV[cellId];		
			const Foam::scalar rc = Foam::pow( (3.0/4.0/Pi)*vc,0.33333333);
			const Foam::scalar phi = 0.579*Foam::pow(vp/vc, 0.132);
			const Foam::scalar sigma2_p = Foam::pow(phi*rp,2);
			const Foam::scalar sigma2_c = Foam::pow(phi*rc,2);

			const auto mu_c = Foam::mag(cp - cellC[cellId]);
			const auto mu2_c = mu_c*mu_c;

			const Foam::scalar delta = 
			(
				sigma2_p * sigma2_p * mu2_c +
				sigma2_p *
				(
					mu2_c + 2*sigma2_c * Foam::log(Foam::sqrt(sigma2_c/sigma2_p)*vp/vc)
				) * (sigma2_c-sigma2_p)
			);

			if(delta<0.0)
			{
				parWeights.push_back({cellId,vp});
				pSubTotal += vp;
				break;
			}

			auto xmax = (-sigma2_p*mu_c +Foam::sqrt(delta))/(sigma2_c-sigma2_p);
			auto xmin = (-sigma2_p*mu_c -Foam::sqrt(delta))/(sigma2_c-sigma2_p);

			auto vpi = 
				0.5 * vp * 
				(
					2 + 
					Foam::erf(xmin/sqrt(2*sigma2_p))-
					Foam::erf(xmax/sqrt(2*sigma2_p))
				)
				+
				0.5 * vc *
				(
					Foam::erf( (xmax-mu_c)/Foam::sqrt(2*sigma2_c) )-
					Foam::erf( (xmin-mu_c)/Foam::sqrt(2*sigma2_c) )
				);
			parWeights.push_back({cellId,vpi});
			pSubTotal += vpi;
		}

		pSubTotal = Foam::max(pSubTotal, static_cast<Foam::scalar>(1.0e-10));
		for(auto& [cellid, w]:parWeights) w /= pSubTotal;
		
	}
	
}



