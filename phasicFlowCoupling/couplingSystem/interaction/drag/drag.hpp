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

	void setSuSpToZero();

	virtual 
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)=0;

	inline
	Foam::scalar dimlessDragI(Foam::scalar Re, Foam::scalar ep) 
	{

		auto Rec = Foam::max(Re,residualRe_);
		Foam::scalar xi = 3.7 - 0.65*Foam::exp(-0.5*Foam::pow(1.5-Foam::log10(Rec),2));
		Foam::scalar Cd = Foam::pow(0.63+4.8/Foam::sqrt(Rec),2);

		return Cd/24 * Re * Foam::pow(ep, -xi ); 
		
	}

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

	
	void calculateDragForce(
		const Plus::realx3ProcCMField& velocity,
		const Plus::realProcCMField& diameter,
		Plus::realx3ProcCMField& particleForce);



	static
	uniquePtr<drag> create(
		Foam::dictionary 		dict, 
		porosity& 				prsty);
		
	
}; 

} // pFlow::coupling


#endif
