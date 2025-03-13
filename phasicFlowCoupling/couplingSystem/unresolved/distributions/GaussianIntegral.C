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
#include <set>

// from OpenFOAM
#include "GaussianIntegral.hpp"
#include "couplingMesh.hpp"

void pFlow::coupling::GaussianIntegral::parseNeighbors(
	const Foam::label 		targetCelli,
	const Foam::vector& 	targetCellCentre, 
	const Foam::scalar 		searchLen,
	const Foam::scalar 		celli,
	std::set<Foam::label>& 	finalList,
	const Foam::label 		layerNumber)
{
	const Foam::labelList& neigborCells = mesh_.cellCells(celli);
	const Foam::vectorField& cellC = mesh_.cellCentres(); 
	const Foam::scalar cellV = mesh_.cellVolumes()[targetCelli];
	const Foam::scalar lCell = 0.5* Foam::pow(cellV, 0.33333);
	// first find all the cells around 
	int found = 0; 
	for(auto nghbrCelli:neigborCells)
	{
		if(nghbrCelli==targetCelli || finalList.count(nghbrCelli) == 1) continue;
		
		Foam::scalar dist = Foam::mag(cellC[nghbrCelli] - targetCellCentre);	
		if(dist <= searchLen + lCell)
		{
			found++;
			finalList.insert(nghbrCelli);
		}
	}

    if(layerNumber <=4 && found != 0)
    {
    	for(auto nghbrCelli:neigborCells)
    	{
    		parseNeighbors
    		(
    			targetCelli,
    			targetCellCentre,
    			searchLen,
    			nghbrCelli,
    			finalList,
    			layerNumber+1
    		);
    	}
    }
}

void pFlow::coupling::GaussianIntegral::constructLists(const Foam::scalar searchLen)
{
	
	if(!mesh_.hasCellCells())
	{
		mesh_.cellCells();
	}
	const Foam::vectorField& cellC = mesh_.cellCentres();
	const Foam::label nCells = mesh_.nCells();
	
	// check for capacity and size 
	if(nCells != neighborList_.size() )
	{
		neighborList_.clear();
		neighborList_.setSize(nCells);
	}
	
	#pragma omp parallel for
	for (Foam::label celli = 0; celli < nCells; celli++) 
	{

        // Get the face labels of the current cell
        auto targetCellCentre = cellC[celli];
        std::set<Foam::label> finalList;

        finalList.insert(celli);

        parseNeighbors(celli, targetCellCentre, searchLen, celli, finalList, 1);        
         
        neighborList_[celli].setSize(finalList.size());
        Foam::label n = 0;
        for(auto nbr:finalList)
        {
        	neighborList_[celli][n]=nbr;
        	n++;
        }
    }
}


pFlow::coupling::GaussianIntegral::GaussianIntegral
(
	Foam::dictionary 		dict, 
	const couplingMesh& 	cMesh,
	const Plus::centerMassField& centerMass
)
:
	neighborLength_(lookupDict<Foam::scalar>(dict, "neighborLength")),
	lengthExtent_(lookupOrDefaultDict(dict, "lengthExtent", static_cast<Foam::scalar>(3.0))),
	weights_("weights",centerMass),
	mesh_(cMesh.mesh())
{
	constructLists(lengthExtent_*neighborLength_);
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


	for(size_t i=0; i<numPar; i++)
	{
		const Foam::label targetCellId = parCellIndex[i];
		auto& weightsI = weights_[i];
		
		weightsI.clear();
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const Foam::labelList& neighbors = neighborList_[targetCellId];

		const Foam::scalar rp = parDiameter[i]/2;
		const Foam::scalar vp = 4*Pi/3 * Foam::pow(rp,3);
		const realx3& cp_i = centerMass[i];
		const Foam::vector cp{cp_i.x(), cp_i.y(), cp_i.z()};

		Foam::scalar pSubTotal = 0;
		for(auto j:neighbors)
		{

			const auto vc = cellV[j];		
			const Foam::scalar rc = Foam::pow( (3.0/4.0/Pi)*vc,0.33333333);
			const Foam::scalar phi = 0.4; //0.579*Foam::pow(vp/vc, 0.132);
			const Foam::scalar sigma2_p = Foam::pow(phi*rp,2);
			const Foam::scalar sigma2_c = Foam::pow(phi*rc,2);

			const auto mu_c = Foam::mag(cp - cellC[j]);
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
				weightsI.push_back(vp);
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
			weightsI.push_back(vpi);
			pSubTotal += vpi;
		}

		pSubTotal = Foam::max(pSubTotal, static_cast<Foam::scalar>(1.0e-15));
		for(auto& w:weightsI) w /= pSubTotal;
		
	}
	

	
}



