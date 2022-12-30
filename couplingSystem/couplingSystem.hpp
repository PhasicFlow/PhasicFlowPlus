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

// from phasicFlow
#include "box.hpp"

#include "scatteredCommunication.hpp"
#include "procCMFields.hpp"

#include "couplingMesh.hpp"
#include "procDEMSystem.hpp"


namespace pFlow::coupling
{

class couplingSystem
{
protected:

	// this is created only on the master processor
	couplingMesh 				couplingMesh_;
	
	MPI::procCommunication 		processorComm_;

	MPI::scatteredCommunication<real> 	realScatteredComm_;

	MPI::scatteredCommunication<realx3> realx3ScatteredComm_;

	MPI::procDEMSystem 			procDEMSystem_;

	MPI::centerMassField 		centerMass_;

	MPI::realProcCMField		particleDiameter_;

	MPI::realx3ProcCMField 		fluidForce_;

	MPI::realx3ProcCMField   	fluidTorque_;

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

	bool collectFluidForce();

	bool collectFluidTorque();

	bool checkParticleDistribution();

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
	auto& fluidTorque()
	{
		return fluidTorque_;
	}

	inline
	auto& fluidForce()
	{
		return fluidForce_;
	}

}; 

} // pFlow::coupling

#endif //__couplingSystem_hpp__
