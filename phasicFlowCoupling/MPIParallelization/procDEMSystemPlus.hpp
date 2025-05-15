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


namespace pFlow
{
class Timers;
class Timer;
}

namespace pFlow::Plus
{

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
		char* argv[]);

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

};

}



#endif //__procDEMSystemPlus_hpp__
