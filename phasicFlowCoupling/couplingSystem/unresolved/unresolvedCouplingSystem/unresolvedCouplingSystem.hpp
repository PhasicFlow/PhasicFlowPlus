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

#ifndef __unresolvedCouplingSystem_hpp__
#define __unresolvedCouplingSystem_hpp__


#include "couplingSystem.hpp"
#include "virtualConstructor.hpp"


namespace pFlow::coupling
{


class unresolvedCouplingSystem
:
	public couplingSystem
{

public:

	unresolvedCouplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	unresolvedCouplingSystem(const unresolvedCouplingSystem&) = delete;
	
	unresolvedCouplingSystem& operator=(const unresolvedCouplingSystem&) = delete;

	unresolvedCouplingSystem(unresolvedCouplingSystem&&) = delete;

	unresolvedCouplingSystem& operator=(unresolvedCouplingSystem&&) = delete;

	~unresolvedCouplingSystem() override = default;

	create_vCtor
	(
		unresolvedCouplingSystem,
		word,
		(
			word demSystemName, 
			Foam::fvMesh& mesh,
			int argc, 
			char* argv[]
		),
		(demSystemName, mesh, argc, argv)
	);

	const Foam::dictionary& unresolvedDict()const
	{
		return this->subDict("unresolved");
	}

	virtual 
	const Foam::volScalarField& alpha()const =0;

	virtual
	const Foam::volScalarField& Sp()const =0;
	
	virtual
	const Foam::volVectorField& Su()const =0;

	virtual
	void calculateFluidInteraction() =0;

	virtual
	void calculatePorosity() =0;

	virtual 
	word shapeTypeName() const= 0;

	static
	uniquePtr<unresolvedCouplingSystem> create
	(
		word demSystemName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]
	);

}; 

} // pFlow::coupling

#endif //__unresolvedCouplingSystem_hpp__
