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
#include "fvCFD.H"

#include "powerDistributionKernel.hpp"


pFlow::coupling::powerDistributionKernel::powerDistributionKernel(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	statistical(dict, cMesh, centerMass, parDiam),
	filterEmpty_(dict.lookupOrDefault<Foam::Switch>("filterEmpty", false))
{
	cellNeighborsSearch();
}


bool pFlow::coupling::powerDistributionKernel::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		 Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		 	);
	
	auto& solidVol = solidVoldTmp.ref();
	size_t numPar = centerMass_.size();
	const vectorField& allCellCntr = this->mesh().cellCentres();
	
	const auto b2 = Foam::pow(neighborLength(),2);

	// lambda for distribution function 
	auto distFunc  = [b2](const Foam::vector& x) -> Foam::scalar
	{
		auto ksi2 = dot(x,x)/b2;
		if(ksi2<1.0)
			return Foam::pow(1-ksi2,4);
		else
			return 0;
	};

	auto distFunc2 = [b2](const Foam::vector& x, const Foam::scalar dist)-> Foam::scalar
	{
		auto ksi2 = dot(x,x)/b2;
		if(ksi2<1.0)
		{
			auto shifted_ksi2 = Foam::min( (dot(x,x)+4*dist*dist)/b2, 1);
			//Foam::Info<<"ksi2 "<< ksi2 <<" shifted_ksi2 "<< shifted_ksi2<<Foam::endl;
			return Foam::pow(1-ksi2,4)+Foam::pow(1-shifted_ksi2,4);
		}
		else
			return 0;
	};
	
	
	//#pragma omp parallel for 
	for(size_t i=0; i<numPar; i++)
	{
		real parVolume = static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));

		const Foam::label targetCellId = parCellIndex_[i];
		
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const auto& nList = neighborList_[targetCellId];
		
		// center of particle 
		const realx3& cp = centerMass_[i];
		Foam::vector CP{cp.x(), cp.y(), cp.z()};	


		std::vector<Foam::scalar> ks; 
		ks.clear();
		ks.reserve(nList.size());
		
		Foam::scalar pSubTotal = 0;
		Foam::label bndryIndex = boundaryCell_[targetCellId].first;
		Foam::label faceIndex = boundaryCell_[targetCellId].second;
		

		if( bndryIndex == -1 )
		{
			// this is not a boundary cell 
			for(auto j:nList)
			{
				Foam::scalar k = distFunc(CP - allCellCntr[j]);
				ks.push_back(k);
				pSubTotal += k;
				
				//Foam::Info<< targetCellId <<" "<< j<<Foam::endl;
			}


			
		}
		else
		{
			// this is a boundary cell 
			const auto& bndry = cMesh_.mesh().boundary()[bndryIndex];
			const Foam::vector parFace = CP - bndry.Cf()[faceIndex];
			const Foam::vector normal = normalised(bndry.Sf()[faceIndex]);
			Foam::scalar parFaceDist = std::abs(normal & parFace);

			for(auto j:nList)
			{	
				Foam::scalar k;
				if(bndryIndex== boundaryCell_[j].first )
				{
					k = distFunc2(CP - allCellCntr[j], parFaceDist);
					//Foam::Info<<"cell "<< targetCellId <<" and j "<< j<<Foam::endl;
				}
				else
					k = distFunc(CP - allCellCntr[j]);
				ks.push_back(k);
				pSubTotal += k;
			}

		}
		
		size_t n = 0;
		pSubTotal = Foam::max(pSubTotal, Foam::vSmall);
		for(auto j:nList)
		{
			solidVol[j] += ks[n++] /pSubTotal * parVolume;
		}

	}
	
	this->ref() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()));

	return true;
}



