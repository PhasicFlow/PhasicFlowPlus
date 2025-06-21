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

#include "diffusion.hpp"
#include "self.hpp"
#include "unresolvedCouplingSystem.hpp"

Foam::tmp<Foam::fvMatrix<Foam::scalar>> pFlow::coupling::diffusion::fvmDdt
(
    const Foam::volScalarField& sField
)
{
	Foam::tmp<Foam::fvMatrix<Foam::scalar>> tfvm
    (
        new Foam::fvMatrix<Foam::scalar>
        (
            sField,
            sField.dimensions()*Foam::dimVol/Foam::dimTime
        )
    );

    Foam::fvMatrix<Foam::scalar>& fvm = tfvm.ref();

    const Foam::scalar rDeltaT = 1.0/dt_.value();

    fvm.diag() = rDeltaT*sField.mesh().Vsc();

    if (sField.mesh().moving())
    {
        fvm.source() = rDeltaT*sField.oldTime().primitiveField()*sField.mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT*sField.oldTime().primitiveField()*sField.mesh().Vsc();
    }

    return tfvm;
}

pFlow::coupling::diffusion::diffusion(
	const unresolvedCouplingSystem& CS,
	const couplingMesh& 			cMesh,
	const Plus::realProcCMField& 	parDiam)
:
	porosity(CS, cMesh, parDiam),
	nSteps_
	(
		Foam::max(1, lookupDict<Foam::label>(CS.unresolvedDict().subDict("porosity"), "nSteps"))
	),
	boundLength_
	(
		lookupDict<Foam::scalar>(CS.unresolvedDict().subDict("porosity"), "boundLength")
	),
	intTime_(boundLength_*boundLength_/4),
	dt_("dt", Foam::dimTime, intTime_/nSteps_),
	DT_("DT", Foam::dimDynamicViscosity/Foam::dimDensity, 1.0),
	picSolDict_("picSolDict")
{
	
	picSolDict_.add("relTol", 0);
	picSolDict_.add("tolerance", 1.0e-8);
	picSolDict_.add("solver", "smoothSolver");
	picSolDict_.add("smoother", "symGaussSeidel");
}


bool pFlow::coupling::diffusion::internalFieldUpdate()
{

	self selfCellDist;
	
	auto solidVolTmp = calculateSolidVol(selfCellDist);

    Foam::fieldRef(*this) = Foam::max(
        1 - solidVolTmp/this->mesh().V(), 
        static_cast<Foam::scalar>(this->alphaMin()) );

	auto picAlphaTmp = Foam::volScalarField::New(
		"picAlpha",
		this->mesh(),
		Foam::dimensioned("picAlpha", Foam::dimless, Foam::scalar(0)),
		"zeroGradient"
	);

	Foam::volScalarField& picAlpha = picAlphaTmp.ref();
	
	Foam::fieldRef(picAlpha) = Foam::min(solidVolTmp/this->mesh().V(), 1.0);
	picAlpha.correctBoundaryConditions();
	
	
	// start of Time loop
	for(Foam::label i=0; i<nSteps_; i++)
	{
		picAlpha.storeOldTime();
		Foam::fvScalarMatrix alphaEq
		(
			fvmDdt(picAlpha) - Foam::fvm::laplacian(DT_,picAlpha)
		);
		alphaEq.solve(picSolDict_);
	}

    //Info<<"*********************************************************\n";
	
	Foam::fieldRef(*this) = 1.0-picAlpha.internalField();

	return true;
}
