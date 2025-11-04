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

#ifndef __solidAveraging_hpp__
#define __solidAveraging_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"


// from phasicFlowPlus
#include "procCMFields.hpp"
#include "couplingMesh.hpp"
#include "porosity.hpp"
#include "schedule.hpp"


namespace pFlow::coupling
{

template<typename DistributionType>
class solidAveraging
{
protected:

    /// averaged property on particle center
    Plus::realx3ProcCMField&  parAvFieldRef_;

public:

    // type info
	TypeInfoTemplate11("solidAveraging", DistributionType);

    solidAveraging(
        const word& type,
        const word& distType,
        const couplingMesh& cMesh,
        Plus::realx3ProcCMField& particleField);

    virtual 
    ~solidAveraging()=default;

    /// virtual constructor  
    create_vCtor
	(
		solidAveraging,
		word,
		(
			const word& type,
            const word& distType,
            const couplingMesh& cMesh,
            Plus::realx3ProcCMField& particleField
		),
		(type, distType, cMesh, particleField)
	);


    /// @brief perform mass-based averaging (if applicable)
    /// @param distributor the cellDistribution system
    /// @param prsty porosity object 
    virtual 
    void calculate(
        const Plus::realx3ProcCMField& particleField
        const DistributorType& cellDistributor,
        const porosity& prsty) = 0;
    
    /// @brief perform number-based averaging (if applicable)
    /// @param distributor the cellDistribution system
    /// @param prsty porosity object 
    virtual 
    void calculateNumberBased(
        const Plus::realx3ProcCMField& particleField
        const DistributorType& distributor
        const porosity& prsty) = 0;
    
    virtual 
    bool requireCellDistribution()const = 0;

    inline
    const realx3& value(Foam::scalar celli, size_t parIdx)const
    {
        return parAvFieldRef_[parIdx];
    }

    inline
    const Plus::realx3ProcCMField& averagedField()const
    {
        return parAvFieldRef_;
    }

    uniquePtr<solidAveraging> create
    (
        const word& type,
        const word& distType,
        const couplingMesh& cMesh,
        Plus::realx3ProcCMField& particleField
    );

};

}


#inlcude "solidAveraging.hpp"

#endif //__solidAveraging_hpp__


/*{
    
        if(!distribute_) return;

        const auto&  parCellInd = prsty.parCellIndex();
        const size_t numPar = parCellInd.size();        
        const auto&  parDiam =  prsty.particleDiameter();
        const auto&  cellVol = prsty.mesh().V();
        const auto&  alpha = prsty.alpha();
        
        
        forAll(Us_,celli)
        {
            Us_[celli] = Foam::Zero;
        }

        #pragma omp parallel for schedule (dynamic)
        for(size_t i=0; i<numPar; i++)
        {
            Foam::scalar pVol = pFlow::Pi/6 *
                    Foam::pow(parDiam[i], static_cast<real>(3.0));

            const Foam::label cellId = parCellInd[i];
            if( cellId >= 0 )
            {
                Foam::vector upv{ pVol*Up_[i].x(), pVol*Up_[i].y(), pVol*Up_[i].z()};
                distributor.distributeValue_OMP(i, cellId, Us_, upv);				
            }
        }

        forAll(Us_,celli)
        {
            Us_[celli] /= Foam::max( (1-alpha[celli])*cellVol[celli], Foam::SMALL);
        }

        distributor.smoothenField(Us_);
        Us_.correctBoundaryConditions();

    }*/