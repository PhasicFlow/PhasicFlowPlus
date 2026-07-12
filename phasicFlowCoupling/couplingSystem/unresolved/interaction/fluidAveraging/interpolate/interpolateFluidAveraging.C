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

#include "interpolateFluidAveraging.hpp"
#include "unresolvedCouplingSystem.hpp"


pFlow::coupling::interpolateFluidAveraging::interpolateFluidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const word&                     name
)
:
    fluidAveraging(type, uCS, name)
{
}

void pFlow::coupling::interpolateFluidAveraging::calculate(const Foam::volVectorField &orgField)
{
    const size_t numPar = averagedField_.size();
    const Plus::centerMassField& centerMass = averagedField_.centerMass();
    const Foam::fvMesh& mesh = uCS_.cMesh().mesh();
    const Plus::procCMField<Foam::label>& parCellIndex = uCS_.parCellIndex();
        
    if(!mesh.hasCellCentres())
    {
        mesh.cellCentres();
    }

    if(!mesh.hasCellCells())
    {
        mesh.cellCells();
    }
    
    const auto& cellCentres = mesh.cellCentres();
    
    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        
        const auto celli = parCellIndex[i];
        if(celli == -1)
        {
            averagedField_[i] = Foam::vector(0,0,0);
            continue;
        } 

        Foam::vector p{centerMass[i].x(), centerMass[i].y(), centerMass[i].z()};
        const auto& cellCentre = cellCentres[celli];

        Foam::scalar dist = Foam::max(p.dist(cellCentre), SMALL);
        const Foam::labelList& neigborCells = mesh.cellCells(celli);
        
        Foam::scalar sum = 1.0/dist;
        Foam::vector sumVel = orgField[celli] * sum;

        forAll(neigborCells, j)
        {
            const auto cellj = neigborCells[j];
            const auto& nCellCentre = cellCentres[cellj];
            
            Foam::scalar nDist = Foam::max(p.dist(nCellCentre), SMALL);
            Foam::scalar w = 1.0/nDist;
            
            sum += w;
            sumVel += orgField[cellj] * w;
        }

        averagedField_[i] =  sumVel / sum;        
    }
}
