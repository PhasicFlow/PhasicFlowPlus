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

#ifndef __drag_hpp__ 
#define __drag_hpp__

// from OpenFOAM
#include "dictionary.H"

// from phasicFlow
#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "porosity.hpp"

namespace pFlow::coupling
{
class drag

{
protected:

	real 					residualRe_;

	porosity& 				porosity_;

	const Foam::volScalarField& 	p_;

	const Foam::volVectorField&		U_;

	Foam::volVectorField 	Su_;

	Foam::volScalarField 	Sp_;

	bool 					isCompressible_ = false;

public:

	// type info
	TypeInfo("drag");

	drag(
		Foam::dictionary 		dict, 
		porosity& 				prsty);

	virtual ~drag() = default;

	create_vCtor
	(
		drag,
		dictionary,
		(
			Foam::dictionary 		dict, 
			porosity& 				prsty
		),
		(dict, prsty)
	);

	Foam::tmp<Foam::volVectorField> 
	pressureGradient(const Foam::volScalarField& rho)const;

	const auto& Su()const
	{
		return Su_;
	}

	const auto& Sp()const
	{
		return Sp_;
	}
	
	inline
	bool isCompressible()const
	{
		return isCompressible_;
	}

	virtual
	void calculateDragForce(
		const MPI::realx3ProcCMField& velocity,
		const MPI::realProcCMField& diameter,
		MPI::realx3ProcCMField& particleForce) = 0;



	static
	uniquePtr<drag> create(
		Foam::dictionary 		dict, 
		porosity& 				prsty);
		
	
}; 

} // pFlow::coupling


#endif
