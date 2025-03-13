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
#ifndef __self_hpp__ 
#define __self_hpp__


// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "typeInfo.hpp"

// from phasicFlowPlus
#include "procCMField.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class self
{
public:

	/// Type info
	TypeInfoNV("self");

	// default constructor 
	self(){}

	/// Construc from dictionary 
	self(
		Foam::dictionary 		dict, 
		const couplingMesh& 	cMesh,
		const Plus::centerMassField& centerMass);

	/// Destructor
	~self() = default;

	
	void updateWeights(
		const Plus::procCMField<Foam::label> & parCellIndex,
		const Plus::procCMField<real> & parDiameter);

	inline
	void distributeValue_OMP(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volScalarField::Internal& internalField,
		const Foam::scalar& val )const
	{
		
		#pragma omp atomic
		internalField[parCellIndx] += val;
		
	}

	inline
	void distributeValue_OMP(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volVectorField::Internal& internalField,
		const Foam::vector& val )const
	{
		
		auto& tv = internalField[parCellIndx]; 
		#pragma omp atomic
		tv.x() += val.x();
		
		#pragma omp atomic
		tv.y() += val.y();

		#pragma omp atomic
		tv.z() += val.z();
		
	}

	inline 
	void distributeValue(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volScalarField::Internal& internalField,
		const Foam::scalar& val )const
	{
		
		internalField[parCellIndx] += val;
		
	}

	inline
	void distributeValue(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volVectorField::Internal& internalField,
		const Foam::vector& val )const
	{
		internalField[parCellIndx] += val;
	}

}; 

} // pFlow::coupling


#endif
