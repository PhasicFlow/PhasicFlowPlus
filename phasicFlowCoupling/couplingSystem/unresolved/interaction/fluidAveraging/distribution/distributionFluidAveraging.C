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

#include "distributionFluidAveraging.hpp"
#include "unresolvedCouplingSystem.hpp"


pFlow::coupling::distributionFluidAveraging::distributionFluidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const word&                     name
)
:
    fluidAveraging(type, uCS, name)
{
}

void pFlow::coupling::distributionFluidAveraging::calculate(const Foam::volVectorField &orgField)
{
    const size_t numPar = averagedField_.size();
    const auto& parCellIndex = uCS_.parCellIndex();
    const auto& distributor = uCS_.distribution();

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        if(parCellIndex[i] != -1)
        {
            distributor.inverseDistributeValue
            (
                i,
                parCellIndex[i],
                orgField,
                averagedField_[i]
            );
        }
    }
}
