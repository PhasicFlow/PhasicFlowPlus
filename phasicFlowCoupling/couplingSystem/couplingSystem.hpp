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

#ifndef __couplingSystem_hpp__
#define __couplingSystem_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "box.hpp"

// from phasicFlow-coupling
#include "particleMapping.hpp"
#include "procCMFields.hpp"
#include "procDEMSystemPlus.hpp"
#include "couplingMesh.hpp"
#include "Timers.hpp"

namespace pFlow::coupling
{


class couplingSystem
:
	public Foam::IOdictionary
{
private:
	
	particleMapping 			particleMapping_;

	couplingMesh 				couplingMesh_;

	Plus::procDEMSystem 		procDEMSystem_;

	Timers 						couplingTimers_;

	Timers 						cfdTimers_;

	Timer 						getDataTimer_;

	Timer 						sendDataTimer_;

	Plus::uint32ProcCMField     particleID_;

	Plus::realProcCMField		particleDiameter_;

	Plus::realx3ProcCMField 	particleVelocity_;

	Plus::realx3ProcCMField 	particleRVelocity_;

	Plus::realx3ProcCMField 	fluidForce_;

	Plus::realx3ProcCMField   	fluidTorque_;

	bool requireRVel_;

	bool collectFluidForce();

	bool collectFluidTorque();

	bool distributeParticles();

protected:

	virtual
	bool distributeParticleFields();

	inline
	Plus::procDEMSystem& pDEMSystem()
	{
		return procDEMSystem_;
	}

	inline 
	particleMapping& parMapping()
	{
		return particleMapping_;
	}

	inline
	Timer& getDataTimer()
	{
		return getDataTimer_;
	}

	inline 
	Timer& sendDataTimer()
	{
		return sendDataTimer_;
	}

	inline 
	Plus::realx3ProcCMField& fluidTorque()
	{
		return fluidTorque_;
	}

	inline
	Plus::realx3ProcCMField& fluidForce()
	{
		return fluidForce_;
	}

	inline
	const Foam::dictionary& dict()const
	{
		return static_cast<const dictionary&>(*this);
	}

public:

	TypeInfo("couplingSystem");

	couplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[],
		bool requireRVel = false);

	couplingSystem(const couplingSystem&) = delete;
	
	couplingSystem& operator=(const couplingSystem&) = delete;

	couplingSystem(couplingSystem&&) = delete;

	couplingSystem& operator=(couplingSystem&&) = delete;

	virtual 
	~couplingSystem() = default;

	virtual
	bool getDataFromDEM(real t, real dt);

	virtual
	bool sendDataToDEM(real t, real dt);

	void sendFluidForceToDEM();

	void sendFluidTorqueToDEM();
	
	bool iterate(real upToTime, bool writeTime, const word& timeName);
	
	inline
	auto& cMesh()
	{
		return couplingMesh_;
	}
	
	inline
	const auto& cMesh()const
	{
		return couplingMesh_;
	}

	inline 
	const particleMapping& parMapping()const
	{
		return particleMapping_;
	}

	inline 
	auto numParticles()const
	{
		return centerMass().size();
	}

	inline 
	auto& particleDiameter()
	{
		return particleDiameter_;
	}

	inline
	Plus::centerMassField& centerMass()
	{
		return particleMapping_.centerMass();
	}

	inline
	const Plus::centerMassField& centerMass()const
	{
		return particleMapping_.centerMass();
	}

	inline 
	Plus::uint32ProcCMField& particleID()
	{
		return particleID_;
	}

	inline
	Plus::realx3ProcCMField& particleVelocity()
	{
		return particleVelocity_;
	}

	inline
	Plus::realx3ProcCMField& particleRVelocity()
	{
		return particleRVelocity_;
	}

	inline 
	const Plus::procVector<box>& meshBoxes()const
	{
		return particleMapping_.meshBoxes();
	}

	inline 
	Timers& cfdTimers()
	{
		return cfdTimers_;
	}

	inline
	Timers& couplingTimers()
	{
		return couplingTimers_;
	}

}; 

} // pFlow::coupling

#endif //__couplingSystem_hpp__
