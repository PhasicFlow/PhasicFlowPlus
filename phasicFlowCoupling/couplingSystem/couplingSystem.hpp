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
#include "scatteredCommunicationPlus.hpp"
#include "procCMFields.hpp"
#include "procDEMSystemPlus.hpp"
#include "couplingMesh.hpp"
#include "Timers.hpp"

namespace pFlow::coupling
{


class couplingSystem
:
	public Foam::IOdictionary,
	public Plus::procCommunication
{
private:
	
	couplingMesh 							couplingMesh_;

	Plus::procVector<box> 					meshBoxes_;

	Plus::scatteredCommunication<real> 		realScatteredComm_;

	Plus::scatteredCommunication<realx3> 	realx3ScatteredComm_;

	Plus::procDEMSystem 		procDEMSystem_;

	Timers 						couplingTimers_;

	Timers 						cfdTimers_;

	Timer 						getDataTimer_;

	Timer 						sendDataTimer_;

	Plus::centerMassField 		centerMass_;

	Plus::uint32ProcCMField     particleID_;

	Plus::realProcCMField		particleDiameter_;

	Plus::realx3ProcCMField 	particleVelocity_;

	Plus::realx3ProcCMField 	particleRVelocity_;

	Plus::realx3ProcCMField 	fluidForce_;

	Plus::realx3ProcCMField   	fluidTorque_;

	bool collectFluidForce();

	bool collectFluidTorque();

	bool distributeParticles();

protected:

	virtual
	bool distributeParticleFields();

	inline
	auto& realx3MappedComm() 
	{
		return realx3ScatteredComm_;
	}

	inline 
	auto& realMappedComm()
	{
		return realScatteredComm_;
	}

	inline
	Plus::procDEMSystem& pDEMSystem()
	{
		return procDEMSystem_;
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
		char* argv[]);

	couplingSystem(const couplingSystem&) = delete;
	
	couplingSystem& operator=(const couplingSystem&) = delete;

	couplingSystem(couplingSystem&&) = delete;

	couplingSystem& operator=(couplingSystem&&) = delete;

	virtual 
	~couplingSystem() = default;

	bool updateMeshBoxes();

	virtual
	bool getDataFromDEM(real t, real dt);

	virtual
	bool sendDataToDEM(real t, real dt);

	/*virtual
	void calculateFluidInteraction() =0;*/

	void sendFluidForceToDEM();

	void sendFluidTorqueToDEM();
	
	bool iterate(real upToTime, bool writeTime, const word& timeName);
	
	inline
	auto& cMesh()
	{
		return couplingMesh_;
	}	

	inline 
	auto numParticles()const
	{
		return centerMass_.size();
	}

	inline 
	auto& particleDiameter()
	{
		return particleDiameter_;
	}

	inline
	Plus::centerMassField& centerMass()
	{
		return centerMass_;
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
		return meshBoxes_;
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
