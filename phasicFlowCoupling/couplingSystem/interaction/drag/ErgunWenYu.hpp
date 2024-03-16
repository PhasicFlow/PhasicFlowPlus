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

#ifndef __ErgunWenYu_hpp__ 
#define __ErgunWenYu_hpp__

#include "drag.hpp"

namespace pFlow::coupling
{
class ErgunWenYu
:
	public drag
{
protected:

	 
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)override
	{
		if( ep >= 0.8 )
		{
			Foam::scalar Cd;
			if(Re <= 1000.0 )
				Cd = 24 * ( 1+0.15*Foam::pow(Re,0.687) ) / Re; 
			else
				Cd = 0.44;
			return Cd/24 * Re * Foam::pow(ep, -3.65 ); 
		}else
		{
			return 	(150.0*(1.0-ep)+ 1.75*Re )/(18.0*ep*ep);
		}

	}

public:

	// type info
	TypeInfo("ErgunWenYu");

	ErgunWenYu(
		Foam::dictionary 		dict, 
		porosity& 				prsty);

	virtual ~ErgunWenYu() = default;

	add_vCtor
	(
		drag,
		ErgunWenYu,
		dictionary
	);

	/*void calculateDragForce(
		const MPI::realx3ProcCMField& velocity,
		const MPI::realProcCMField& diameter,
		MPI::realx3ProcCMField& particleForce)override;*/
		
	
}; 

} // pFlow::coupling


#endif
