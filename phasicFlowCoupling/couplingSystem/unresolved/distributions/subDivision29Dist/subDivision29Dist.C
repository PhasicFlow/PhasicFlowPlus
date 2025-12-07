#include "subDivision29Dist.hpp"
#include "couplingMesh.hpp"

const pFlow::real sin_45[] = {0.7071067811865475,  0.7071067811865475, -0.7071067811865475, -0.7071067811865475};
const pFlow::real cos_45[] = {0.7071067811865475, -0.7071067811865475, -0.7071067811865475,  0.7071067811865475};

pFlow::coupling::subDivision29Dist::subDivision29Dist
(
    const Foam::dictionary& 	 parrentDict, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    distribution(parrentDict, cMesh, centerMass)
{
    
}

void pFlow::coupling::subDivision29Dist::updateWeights(const Plus::procCMField<real> &parDiameter)
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

        bool fullInside;
        bool halfInside;

        const Foam::point 	pPos{centerMass[i].x(),centerMass[i].y(),centerMass[i].z()};
        const Foam::scalar 	pRad = parDiameter[i]/2;
        
                
        cmesh.pointSphereInCell(
            pPos, 
            0.62392*pRad, 
            0.917896*pRad, 
            cntrCellId, 
            halfInside, 
            fullInside);
        
        if(fullInside)
        {
            parWeights.push_back({cntrCellId,1.0});
            continue;
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
            cmesh.findPointsInCells(hpoints, cntrCellId,nCellIds, hcellIds );
            
            parWeights.push_back({cntrCellId,(29-nCellIds)/29.0});
            for(auto ci=0; ci<nCellIds; ci++ )
            {
                parWeights.push_back({hcellIds[ci],1.0/29.0});	
            }

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
            cmesh.findPointsInCells(points, cntrCellId,nCellIds, cellIds);
            
            parWeights.push_back({cntrCellId,(29-nCellIds)/29.0});
            for(auto ci=0; ci<nCellIds; ci++ )
            {
                parWeights.push_back({cellIds[ci],1.0/29.0});	
            }
        }
    }

}
