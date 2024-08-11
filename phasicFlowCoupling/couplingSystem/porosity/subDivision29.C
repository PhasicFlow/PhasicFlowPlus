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

#include "omp.h"

#include "subDivision29.hpp"
#include "streams.hpp"

const pFlow::real sin_45[] = {0.7071067811865475,  0.7071067811865475, -0.7071067811865475, -0.7071067811865475};
const pFlow::real cos_45[] = {0.7071067811865475, -0.7071067811865475, -0.7071067811865475,  0.7071067811865475};



pFlow::coupling::subDivision29::subDivision29(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam)
{

}


bool pFlow::coupling::subDivision29::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		 Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		 	);
	
	auto& solidVol = solidVoldTmp.ref();
	
	size_t numPar = centerMass_.size();

	#pragma omp parallel for 
	for(size_t i=0; i<numPar; i++)
	{

		const Foam::label cntrCellId = parCellIndex_[i];
		if( cntrCellId < 0 )continue;

		Foam::FixedList<realx3, 28> points;
		Foam::FixedList<Foam::label, 28> cellIds;

		const realx3 pPos = centerMass_[i];
		const real 	 pRad = particleDiameter_[i]/2;
		
		// 4*Pi/3 = 4.1887902047864
		const real pSubVol = static_cast<real>(4.1887902047864/29.0) *
					pFlow::pow(pRad, static_cast<real>(3.0));

		realx3 offset(0,0,0);	
		Foam::label n = 0;

		for (real r=0.62392*pRad; r<pRad; r+=static_cast<real>(0.293976*pRad) )
		{
			// 8 subdivisions of particle
			for(int32 i_alp =0; i_alp<4; i_alp++)
			{
				for(int32 i_bet=0; i_bet<2;i_bet++)
				{
					offset = {  
						r*sin_45[i_alp]*cos_45[i_bet],
						r*sin_45[i_alp]*sin_45[i_bet],
						r*cos_45[i_alp] };

					points[n++] = pPos + offset;
				}
			}

			for( int j=-1; j<=1; j+=2 )
			{
	        	offset= {r*j, 0.0, 0.0};
	          	
	          	points[n++] = pPos + offset;
	            
	            offset = {0.0, r*j, 0.0};
	            points[n++] = pPos + offset;
	            

	            offset = {0.0, 0.0, r*j};
	            points[n++] = pPos + offset;
			}
		}

		Foam::label nCellIds = 0;
		cMesh_.findPointsInCells(points, cntrCellId,nCellIds, cellIds );
		
		for(auto ci=0; ci<nCellIds; ci++ )
		{
			#pragma omp atomic
			solidVol[cellIds[ci]] += pSubVol;	
		}

		#pragma omp atomic
		solidVol[cntrCellId] += (29-nCellIds)*pSubVol;
	
	}// omp parallel 
	

	this->ref() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()) );

	return true;
}
