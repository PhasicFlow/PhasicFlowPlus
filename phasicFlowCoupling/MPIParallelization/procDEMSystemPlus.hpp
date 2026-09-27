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

#ifndef __procDEMSystemPlus_hpp__ 
#define __procDEMSystemPlus_hpp__


// from phasicFlow
#include "DEMSystem.hpp"
#include "box.hpp"


// from coupling-phasicFlow
#include "procVectorPlus.hpp"
#include <string>


namespace pFlow
{
class Timers;
class Timer;
}

namespace pFlow::Plus
{

/**
 * @class procDEMSystem
 * @brief Manages top-level DEM execution and parallel MPI memory routing.
 *
 * Serves as the principal interface between OpenFOAM processors and the
 * underlying Kokkos GPU/CPU particle system data streams.
 */
class procDEMSystem
{
protected:

	// this will be nullptr for all processors except 
	// the main processor 
	uniquePtr<DEMSystem> demSystem_ = nullptr;

	real 				 startTime_= 0;
public:

	procDEMSystem(
		word demSystemName,
		int argc, 
		char* argv[],
		bool requireRVel = false);

	virtual ~procDEMSystem() = default;

	// --- global time & control ---

	inline 
	real startTime()const
	{
		return startTime_;
	}

	inline 
	bool getDataFromDEM()
	{
		if(demSystem_)
		{
			return demSystem_->beforeIteration();
		}
		else
		{
			return true;
		}
	}

	
	bool updateParticleDistribution(real extentFraction, const procVector<box>& domains)
	{
		if(demSystem_)
		{
			/*demSystem_->updateParticleDistribution(extentFraction, domains);
			output<<"numProcessos "<< domains.size()<<endl;
			output<<" num in Proc0 "<<demSystem_->numParInDomain(0)<<endl;*/
			
			return demSystem_->updateParticleDistribution(extentFraction, domains);
		}
		else
		{
			return true;
		}
	}

	// --- mechanical field accessors ---

	inline 
	span<const int32> parIndexInDomain(int32 di)const
	{
		if(demSystem_)
		{
			return demSystem_->parIndexInDomain(di);
		}else
		{
			return span<const int32>();
		}
	}

	inline
	span<uint32> particleIdAllMaster()const
	{
		if(demSystem_)
		{
			return demSystem_->particleId();
		}else
		{
			return span<uint32>();
		}
	}

	inline 
	span<realx3> particlesCenterMassAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->position();
		}else
		{
			return span<realx3>();
		}
	}

	inline 
	span<realx3> particlesVelocityAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->velocity();
		}else
		{
			return span<realx3>();
		}
	}

	inline 
	span<realx3> particlesRVelocityAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->rVelocity();
		}else
		{
			return span<realx3>();
		}
	}


	inline 
	span<realx3> particlesFluidForceAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidForce();
		}else
		{
			return span<realx3>();
		}
	}

	inline
	span<realx3> particlesAccelerationAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->acceleration();
		}else
		{
			return span<realx3>();
		}
	}

	inline
	span<realx3> particlesFluidTorqueAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidTorque();
		}else
		{
			return span<realx3>();
		}
	}

	inline
	std::vector<real> shapeDiametersAllMaster()const
	{
		if(demSystem_)
		{
			return demSystem_->shapeDiameters();
		}
		else
		{
			return std::vector<real>{};
		}
	}

	inline
	span<real> particlesDiameterAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->diameter();
		}else
		{
			return span<real>();
		}		
	}

	inline
	span<real> particlesCourseGrainFactorMasterAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->courseGrainFactor();
		}else
		{
			return span<real>();
		}		
	}

	inline
	procVector<int32> numParInDomainMaster()const
	{
		if(demSystem_)
		{
			return demSystem_->numParInDomains();
		}
		else
		{
			return procVector<int32>(true);
		}
	}


	inline
	procVector<span<const int32>> parIndexInDomainsMaster()const
	{
		procVector<span<const int32>> parIndex(true);
		if(demSystem_)
		{
			for(size_t i=0; i<parIndex.size(); i++)
			{
				parIndex[i] = demSystem_->parIndexInDomain(i);
			}
		}

		return parIndex;
	}

	// --- thermal coupling extensions ---

	inline
	span<real> emissivity()
	{
		if(demSystem_)
		{
			return demSystem_->emissivity();
		}else
		{
			return span<real>();
		}
	}

	// True only when a thermal interaction exists AND radiation is
	// enabled in the interaction dictionary. Lets callers skip
	// radSumTemp()/radNumPrt() entirely when radiation is off, rather
	// than distributing a zero-filled buffer nobody uses.
	inline
	bool hasRadiation()const
	{
		if(demSystem_)
		{
			return demSystem_->hasRadiation();
		}else
		{
			return false;
		}
	}

	inline
	span<real> radSumTemp()
	{
		if(demSystem_)
		{
			return demSystem_->radSumTemp();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<uint32> radNumPrt()
	{
		if(demSystem_)
		{
			return demSystem_->radNumPrt();
		}else
		{
			return span<uint32>();
		}
	}

	inline
	span<real> parFluidHeatSourceConv()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidHeatSourceConv();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> parFluidHeatSourceRad()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidHeatSourceRad();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> parFluidKappa()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidKappa();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> parFluidAlpha()
	{
		if(demSystem_)
		{
			return demSystem_->parFluidAlpha();
		}else
		{
			return span<real>();
		}
	}

	// --- multi-species chemical reaction extensions ---
	// Reaction coupling (not yet reviewed).
	/*
	inline
	span<real> solidMassFractions()
	{
		if(demSystem_)
		{
			return demSystem_->solidMassFractions();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> gasMassSource()
	{
		if(demSystem_)
		{
			return demSystem_->gasMassSource();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> gasMassSourceSp()
	{
		if(demSystem_)
		{
			return demSystem_->gasMassSourceSp();
		}else
		{
			return span<real>();
		}
	}

	inline
	span<real> gasConcentrations()
	{
		if(demSystem_)
		{
			return demSystem_->gasConcentrations();
		}else
		{
			return span<real>();
		}
	}

	// Per-particle solid-side reaction heat: (1-eta)*Q_rxn [W].
	inline
	span<real> reactionHeat()
	{
		if(demSystem_)
		{
			return demSystem_->reactionHeat();
		}else
		{
			return span<real>();
		}
	}

	// Per-particle fluid-side reaction heat: eta*Q_rxn [W].
	// Zero-filled when eta = 0 (default for surface reactions).
	inline
	span<real> reactionHeatFluid()
	{
		if(demSystem_)
		{
			return demSystem_->reactionHeatFluid();
		}else
		{
			return span<real>();
		}
	}

	// Returns the DEM gas species names as a standard C++ vector. Clean
	// C++ interface, no OpenFOAM dependencies, to keep standalone DEM
	// capability. Used by buildSpeciesMapping().
	inline
	std::vector<std::string> gasSpeciesNames()const
	{
		if(demSystem_)
		{
			return demSystem_->gasSpeciesNames();
		}else
		{
			return std::vector<std::string>();
		}
	}

	// Gas species molar masses [kg/mol] from the DEM kinetics, in the
	// same order as gasSpeciesNames(). Used by buildSpeciesMapping()'s
	// cross-check at CFD startup to catch a stale/mismatched
	// transportProperties/gasMw entry.
	inline
	std::vector<real> gasMolarMasses()const
	{
		if(demSystem_)
		{
			return demSystem_->gasMolarMasses();
		}else
		{
			return std::vector<real>();
		}
	}

	// sendReactionDataToDEM() has no wrapper here: reaction data flows
	// one way from DEM to CFD via reactionDataHostUpdatedSync() inside
	// getDataFromDEM(), and the other way via sendGasConcentrationsToDEM()
	// below - there is no separate CFD-to-DEM reaction-data send this
	// class needs to expose.
	*/

	// --- sync dispatchers: host -> dem device ---

	inline 
	bool sendFluidForceToDEM()
	{
		if(demSystem_)
		{
			return demSystem_->sendFluidForceToDEM();
		}
		else
		{
			return true;
		}
	}

	inline
	bool sendFluidTorqueToDEM()
	{
		if(demSystem_)
		{
			return demSystem_->sendFluidTorqueToDEM();
		}
		else
		{
			return true;
		}
	}

	inline
	bool sendFluidHeatSourcesToDEM()
	{
		if(demSystem_)
		{
			return demSystem_->sendFluidHeatSourcesToDEM();
		}
		else
		{
			return true;
		}
	}

	inline
	bool sendFluidPropertiesToDEM()
	{
		if(demSystem_)
		{
			return demSystem_->sendFluidPropertiesToDEM();
		}
		else
		{
			return true;
		}
	}

	/*
	inline
	bool sendGasConcentrationsToDEM()
	{
		if(demSystem_)
		{
			return demSystem_->sendGasConcentrationsToDEM();
		}
		else
		{
			return true;
		}
	}
	*/

	// --- timestep execution ---

	inline
	bool iterate(real upToTime, bool writeTime, const word& timeName)
	{
		
		if(demSystem_)
		{
			if(writeTime)
				return demSystem_->iterate(upToTime, upToTime, timeName);
			else
				return demSystem_->iterate(upToTime);
		}
		else
		{
			return true;
		}
	}

	inline
	bool iterate(real upToTime)
	{
		if(demSystem_)
		{
			return demSystem_->iterate(upToTime);
		}
		else
		{
			return true;
		}
	}

	inline
	Timers* getTimers()
	{
		if(demSystem_)
		{
			return &demSystem_->Control().timers();
		}
		else
		{
			return nullptr;
		}
	}

	inline
	span<real> particlesTemperatureAllMaster()
	{
		if(demSystem_)
		{
			return demSystem_->temperature();
		}else
		{
			return span<real>();
		}
	}

};

}



#endif //__procDEMSystemPlus_hpp__
