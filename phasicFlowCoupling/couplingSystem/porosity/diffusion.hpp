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

#ifndef __diffusion_hpp__ 
#define __diffusion_hpp__

#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "PIC.hpp"


namespace pFlow::coupling
{

/**
 * Particle In Cell (diffusion) model for porosity calculation
 * 
 * This model only considers the particle center and if the particle center 
 * resides inside a cell, it is assumed that the whole volume of the particle
 * is located in that cell.
 * 
 */
class diffusion
: 
	public PIC
{
private:

	Foam::label 			  nSteps_;

	Foam::scalar 			  boundLength_;

	Foam::scalar 			  intTime_;

	Foam::dimensionedScalar   dt_;

	Foam::dimensionedScalar   DT_;

	Foam::dictionary          picSolDict_;



	Foam::tmp<Foam::fvMatrix<Foam::scalar>> fvmDdt
    (
        const Foam::volScalarField& sField
    );

public:

	/// Type info
	TypeInfo("diffusion");

	/// Construc from dictionary 
	diffusion(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		Plus::centerMassField& 	centerMass, 
		Plus::realProcCMField& 	parDiam);

	/// Destructor
	virtual ~diffusion() = default;

	/// Add this constructor to the list of virtual constructors
	add_vCtor
	(
		porosity,
		diffusion,
		dictionary
	);

	bool internalFieldUpdate() override;

	
}; 

} // pFlow::coupling


#endif
