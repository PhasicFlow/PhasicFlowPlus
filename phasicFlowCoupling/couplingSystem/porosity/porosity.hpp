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

#ifndef __porosity_hpp__ 
#define __porosity_hpp__

// from OpenFOAM
#include "dictionary.H"
#include "volFields.H"

// from phasicFlow
#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "couplingMesh.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{


class porosity
:
	public Foam::volScalarField
{
protected:

	real 							alphaMin_;

	couplingMesh& 					cMesh_;

	MPI::centerMassField& 			centerMass_;

	MPI::realProcCMField&  			particleDiameter_;

	MPI::procCMField<Foam::label> 	parCellIndex_;

public:

	// type info
	TypeInfo("porosity");

	porosity(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	virtual ~porosity() = default;

	create_vCtor
	(
		porosity,
		dictionary,
		(
			Foam::dictionary		dict, 
			couplingMesh& 			cMesh, 
			MPI::centerMassField& 	centerMass, 
			MPI::realProcCMField& 	parDiam
		),
		(dict, cMesh, centerMass, parDiam)
	);

	
	const Foam::volScalarField& alpha()const
	{
		return *this;
	}

	inline
	real alphaMin()const
	{
		return alphaMin_;
	}

	inline 
	const auto& particleCellIndex()const
	{
		return parCellIndex_;
	}

	inline
	const auto& centerMass()const
	{
		return centerMass_;
	}

	inline 
	const auto& particleDiameter()const
	{
		return particleDiameter_;
	}

	
	virtual
	bool calculatePorosity() = 0;



	static
	uniquePtr<porosity> create(
		Foam::dictionary		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);
	

}; 

} // pFlow::coupling


#endif
