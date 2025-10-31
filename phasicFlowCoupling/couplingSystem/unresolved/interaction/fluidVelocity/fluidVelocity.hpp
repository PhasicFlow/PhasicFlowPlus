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

#ifndef __fluidVelocity_hpp__
#define __fluidVelocity_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"
#include "procCMFields.hpp"
#include "couplingMesh.hpp"

namespace pFlow::coupling
{

//class couplingMesh;

class fluidVelocity
{
    /// fluid velocity at cell center
    const Foam::volVectorField&     Uf_;

    const bool                      interpolate_;

    // interpolat using cell distribution? 
    const bool                      cellDistribution_;

    /// interpolated fluid velocity at particle center of mass
    pFlow::Plus::procCMField<Foam::vector>   Up_;

    void interpolateCells(const couplingMesh& cMesh);

    template<typename DistributorType>
    void interpolateUsingDistributor
    (
        const DistributorType& distributor,
        const couplingMesh& cMesh
    )
    {
        const size_t numPar = Up_.size();
        const auto& parCellIndex = cMesh.parCellIndex();

        #pragma omp parallel for schedule (dynamic)
        for(size_t i=0; i<numPar; i++)
        {
            distributor.inverseDistributeValue
            (
                i,
                parCellIndex[i],
                Uf_,
                Up_[i]
            );
        }
    }

public:

    fluidVelocity(
        const word& type, 
        const Foam::volVectorField& U,
        const couplingMesh& cMesh);

    template<typename DistributorType>
    void interpolate(
        const DistributorType& distributor,
        const couplingMesh& cMesh
    )
    {
        if(!interpolate_)
        {
            return;
        }

        if(cellDistribution_)
        {
            interpolateUsingDistributor(distributor, cMesh);
        }
        else
        {
            interpolateCells(cMesh);
        }
        return ;
    }

    inline
    const Foam::vector& uFluid(Foam::scalar celli, size_t parIdx)const
    {
        if(interpolate_)
        { 
            return Up_[parIdx];
        }
        else
        {
            return Uf_[celli];
        }
    }

    inline 
    bool requireCellDistribution()const
    {
        return cellDistribution_;
    }

};

}


#endif //__fluidVelocity_hpp__
