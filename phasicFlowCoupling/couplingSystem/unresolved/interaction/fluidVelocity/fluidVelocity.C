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


#include "fluidVelocity.hpp"
#include "couplingMesh.hpp"
#include "streams.hpp"

pFlow::coupling::fluidVelocity::fluidVelocity
(
    const word &type, 
    const Foam::volVectorField &U, 
    const couplingMesh &cMesh
)
:
    Uf_(U),
    interpolate_(type == "particle"),
    Up_
    (
        "Up",
        cMesh.centerMass()
    )
{
}

void pFlow::coupling::fluidVelocity::interpolate(const couplingMesh &cMesh)
{
    if(!interpolate_)
    {
        return;
    }
    const auto& mesh = cMesh.mesh();
    if(!mesh.hasCellCentres())
    {
        mesh.cellCentres();
    }
    const auto& cellCentres = mesh.cellCentres();
    const auto& parCellIndex = cMesh.parCellIndex();
    const auto& centerMass = cMesh.centerMass();
    const size_t nPar = centerMass.size();
    
    #pragma ParallelRegion
    for(size_t i=0; i<nPar; i++)
    {
        
        const auto celli = parCellIndex[i];
        if(celli == -1)
        {
            Up_[i] = zero3;
            continue;
        } 

        Foam::vector p{centerMass[i].x(), centerMass[i].y(), centerMass[i].z()};
        const auto& cellCentre = cellCentres[celli];

        Foam::scalar dist = Foam::max(p.dist(cellCentre), SMALL);
        const Foam::labelList& neigborCells = mesh.cellCells(celli);
        
        Foam::scalar sum = 1.0/dist;
        Foam::vector sumVel = Uf_[celli] * sum;

        forAll(neigborCells, j)
        {
            const auto cellj = neigborCells[j];
            const auto& nCellCentre = cellCentres[cellj];
            
            Foam::scalar nDist = Foam::max(p.dist(nCellCentre), SMALL);
            Foam::scalar w = 1.0/nDist;
            
            sum += w;
            sumVel += Uf_[cellj] * w;
        }

        Up_[i] =  realx3(sumVel[0] / sum, 
                         sumVel[1] / sum, 
                         sumVel[2] / sum);
                
    }
    


}
