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
#include "IOdictionary.H"

// from phasicFlow
#include "box.hpp"

// from phasicFlow-coupling
#include "scatteredCommunication.hpp"
#include "procCMFields.hpp"
#include "procDEMSystem.hpp"
#include "couplingMesh.hpp"
#include "porosity.hpp"
#include "drag.hpp"

namespace pFlow::coupling
{


class couplingSystem
:
	public Foam::IOdictionary
{
protected:

	pFlow::real					subDomainExpansionFraction_;

	pFlow::real 				subDomainUpdateInterval_;

	pFlow::real 				lastTimeUpdated_;
	
	couplingMesh 				couplingMesh_;
	
	MPI::procVector<box> 		meshBoxes_;

	MPI::procCommunication 		processorComm_;

	MPI::scatteredCommunication<real> 	realScatteredComm_;

	MPI::scatteredCommunication<realx3> realx3ScatteredComm_;

	MPI::procDEMSystem 			procDEMSystem_;

	MPI::centerMassField 		centerMass_;

	MPI::realProcCMField		particleDiameter_;

	MPI::realx3ProcCMField 		particleVelocity_;

	MPI::realx3ProcCMField 		fluidForce_;

	MPI::realx3ProcCMField   	fluidTorque_;

	uniquePtr<porosity>			porosity_ = nullptr;

	uniquePtr<drag> 			drag_ = nullptr;

	bool checkForDomainUpdate(real t, real fluidDt);

public:


	couplingSystem(
		word demSystemName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	couplingSystem(const couplingSystem&) = delete;
	
	couplingSystem& operator=(const couplingSystem&) = delete;

	couplingSystem(couplingSystem&&) = delete;

	couplingSystem& operator=(couplingSystem&&) = delete;

	virtual ~couplingSystem() =default;

	bool getDataFromDEM(real t, real dt);

	bool calculatePorosity();

	bool collectFluidForce();

	bool collectFluidTorque();

	bool distributeParticles();

	bool distributeVelocity();

	bool iterate(real upToTime, real writeTime, const word& timeName)
	{
		return procDEMSystem_.iterate(upToTime, writeTime, timeName);
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
	auto& centerMass()
	{
		return centerMass_;
	}

	inline
	auto& particleVelocity()
	{
		return particleVelocity_;
	}

	inline 
	auto& fluidTorque()
	{
		return fluidTorque_;
	}

	inline
	auto& fluidForce()
	{
		return fluidForce_;
	}

	inline 
	const auto& meshBoxes()const
	{
		return meshBoxes_;
	}



}; 

} // pFlow::coupling

#endif //__couplingSystem_hpp__
