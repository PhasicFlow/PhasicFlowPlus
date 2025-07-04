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

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "typeInfo.hpp"

// from phasicFlowPlus
#include "procCMField.hpp"
#include "self.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class diffusion
:
    public self
{
    Foam::label 			    nSteps_;

	Foam::scalar 			    standardDeviation_;

	Foam::scalar 			    intTime_;

	Foam::dimensionedScalar     dt_;

	Foam::dimensionedScalar     DT_;

	Foam::dictionary            smoothSolDict_;

    Foam::tmp<Foam::fvMatrix<Foam::scalar>> fvmDdt
    (
        const Foam::volScalarField& sField
    )const;

    Foam::tmp<Foam::fvMatrix<Foam::vector>> fvmDdt
    (
        const Foam::volVectorField& sField
    )const;

public:

    /// Type info
	TypeInfoNV("diffusion");

	/// Construc from dictionary 
	diffusion(
		Foam::dictionary 		dict, 
		const couplingMesh& 	cMesh,
		const Plus::centerMassField& centerMass);

    ~diffusion()=default;

    void smoothenField(Foam::volScalarField& field)const;

    void smoothenField(Foam::volVectorField& field)const;

    bool requireCellDistribution()const 
	{
		return false;
	}
};

}

#endif //__diffusion_hpp__