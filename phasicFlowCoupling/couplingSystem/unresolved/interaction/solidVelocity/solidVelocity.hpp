
#ifndef __solidVelocity_hpp__
#define __solidVelocity_hpp__

#include "OFCompatibleHeader.hpp"
#include "procCMFields.hpp"
#include "couplingMesh.hpp"
#include "porosity.hpp"
#include "schedule.hpp"

namespace pFlow::coupling
{

class solidVelocity
{

    const Plus::realx3ProcCMField&  Up_;

    const bool                      distribute_;

    Foam::volVectorField            Us_;

public:

    solidVelocity(
        const word& type,
        const Plus::realx3ProcCMField& parVel,
        const couplingMesh& cMesh);

    template<typename DistributorType>
    void average(
        const DistributorType& distributor,
        const porosity& prsty)
    {
    
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

        #pragma ParallelRegion
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
        Us_.correctBoundaryConditions();

    }

    inline
    Foam::vector vSolid(Foam::scalar celli, size_t parIdx)const
    {
        if(distribute_)
        {
            return Us_[celli];
        }
        else
        {
            const auto& up = Up_[parIdx]; 
            return Foam::vector{up.x(), up.y(), up.z()};
        }
    }

    inline 
    bool requireCellDistribution()const
    {
        return distribute_;
    }
};

}

#endif //__solidVelocity_hpp__
