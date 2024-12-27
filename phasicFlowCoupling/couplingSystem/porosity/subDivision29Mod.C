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


#include "subDivision29Mod.hpp"
#include "streams.hpp"

const pFlow::real sin_45[] = {0.7071067811865475,  0.7071067811865475, -0.7071067811865475, -0.7071067811865475};
const pFlow::real cos_45[] = {0.7071067811865475, -0.7071067811865475, -0.7071067811865475,  0.7071067811865475};



pFlow::coupling::subDivision29Mod::subDivision29Mod(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	Plus::centerMassField& 	centerMass, 
	Plus::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam)
{

}


bool pFlow::coupling::subDivision29Mod::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		 Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		 	);
	
	auto& solidVol = solidVoldTmp.ref();
		
	size_t numPar = centerMass_.size();
	

	#pragma omp parallel for schedule (dynamic)
	for(size_t i=0; i<numPar; i++)
	{
		const Foam::label cntrCellId = parCellIndex_[i];
		if( cntrCellId < 0 )continue;

		bool fullInside;
		bool halfInside;

		const Foam::point 	pPos{centerMass_[i].x(),centerMass_[i].y(),centerMass_[i].z()};
		const Foam::scalar 	pRad = particleDiameter_[i]/2;
		
		// 4*Pi/3 = 4.1887902047864
		const real pSubVol = static_cast<real>(4.1887902047864/29.0) *
					pFlow::pow(pRad, static_cast<real>(3.0));

		
		cMesh_.pointSphereInCell(
			pPos, 
			0.62392*pRad, 
			0.917896*pRad, 
			cntrCellId, 
			halfInside, 
			fullInside);	

		if(fullInside)
		{
			#pragma omp atomic
			solidVol[cntrCellId] += 29*pSubVol;

		}		
		else if(halfInside)
		{
			Foam::point offset(0,0,0);
			Foam::FixedList<Foam::point, 14> hpoints;
			Foam::FixedList<Foam::label, 14> hcellIds;

			Foam::label n = 0;

			Foam::scalar r = 0.917896*pRad;

			for(int32 i_alp =0; i_alp<4; i_alp++)
			{
				for(int32 i_bet=0; i_bet<2;i_bet++)
				{
					offset = {  
						r*sin_45[i_alp]*cos_45[i_bet],
						r*sin_45[i_alp]*sin_45[i_bet],
						r*cos_45[i_alp] };

					hpoints[n++] = pPos + offset;
				}
			}

			for( int j=-1; j<=1; j+=2 )
			{
	        	offset= {r*j, 0.0, 0.0};
	          	
	          	hpoints[n++] = pPos + offset;
	            
	            offset = {0.0, r*j, 0.0};
	            hpoints[n++] = pPos + offset;
	            

	            offset = {0.0, 0.0, r*j};
	            hpoints[n++] = pPos + offset;
			}


			Foam::label nCellIds = 0;
			cMesh_.findPointsInCells(hpoints, cntrCellId,nCellIds, hcellIds );
			
			for(auto ci=0; ci<nCellIds; ci++ )
			{
				#pragma omp atomic
				solidVol[hcellIds[ci]] += pSubVol;	
			}

			#pragma omp atomic
			solidVol[cntrCellId] += (29-nCellIds)*pSubVol;

			continue;
		}
		else
		{
			Foam::point offset(0,0,0);
			Foam::FixedList<Foam::point, 28> points;
			Foam::FixedList<Foam::label, 28> cellIds;
			Foam::label n=0;
			for (Foam::scalar r=0.62392*pRad; r<pRad; r+=0.293976*pRad) 
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
			cMesh_.findPointsInCells(points, cntrCellId,nCellIds, cellIds);
			
			for(auto ci=0; ci<nCellIds; ci++ )
			{
				#pragma omp atomic
				solidVol[cellIds[ci]] += pSubVol;	
			}

			#pragma omp atomic
			solidVol[cntrCellId] += (29-nCellIds)*pSubVol;
		}
	
	}// omp parallel for
	
	
	this->ref() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()) );

	return true;
}
