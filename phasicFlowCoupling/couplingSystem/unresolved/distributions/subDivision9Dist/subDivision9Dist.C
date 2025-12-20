#include "subDivision9Dist.hpp"
#include "couplingMesh.hpp"

const pFlow::real sin_45[] = {0.7071067811865475,  0.7071067811865475, -0.7071067811865475, -0.7071067811865475};
const pFlow::real cos_45[] = {0.7071067811865475, -0.7071067811865475, -0.7071067811865475,  0.7071067811865475};

pFlow::coupling::subDivision9Dist::subDivision9Dist
(
    const Foam::dictionary& 	 parrentDict, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    distribution(parrentDict, cMesh, centerMass)
{
    
}

void pFlow::coupling::subDivision9Dist::updateWeights(const Plus::procCMField<real> &parDiameter)
{

    const auto& cmesh = this->cMesh();
    const auto& parCellIndex = cmesh.parCellIndex();
    auto& weights = this->weights();
    const auto& centerMass = weights.centerMass();
    const size_t numPar = centerMass.size();

    #pragma omp parallel for schedule (dynamic)
	for(size_t i=0; i<numPar; i++)
	{

		const Foam::label cntrCellId = parCellIndex[i];
        auto& parWeights = weights[i];
        
        parWeights.clear();
        if( cntrCellId < 0 )continue;

		Foam::FixedList<realx3, 8> points;
		Foam::FixedList<Foam::label, 8> cellIds;

		const realx3 pPos = centerMass[i];
		const real pRad = parDiameter[i]/2;
		
		realx3 offset(0,0,0);

		Foam::label n = 0;
		real r = static_cast<real>(0.5*1.48075) * pRad;
		
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

		Foam::label nCellIds = 0;
		cmesh.findPointsInCells(points, cntrCellId, nCellIds, cellIds );
        
        parWeights.push_back({cntrCellId,(9-nCellIds)/9.0});

		for(auto ci=0; ci<nCellIds; ci++ )
		{
            parWeights.push_back({cellIds[ci],1.0/9.0});
		}

	} // omp parallel for

}
