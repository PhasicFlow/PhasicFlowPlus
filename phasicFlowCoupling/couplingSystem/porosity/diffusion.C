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

// from OpenFOAM
#include "fvCFD.H"


#include "diffusion.hpp"


Foam::tmp<Foam::fvMatrix<Foam::scalar>> pFlow::coupling::diffusion::fvmDdt
(
    const Foam::volScalarField& sField
)
{
	Foam::tmp<fvMatrix<Foam::scalar>> tfvm
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
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	PIC(dict, cMesh, centerMass, parDiam),
	nSteps_(Foam::max(1,dict.lookup<Foam::label>("nSteps"))),
	boundLength_(dict.lookup<Foam::scalar>("boundLength")),
	intTime_(0.75*0.75*boundLength_*boundLength_/4),
	dt_("dt", Foam::dimTime, intTime_/nSteps_),
	DT_("DT", Foam::dimDynamicViscosity/Foam::dimDensity, 1.0),
	picSolDict_("picSolDict")
{
	
	picSolDict_.add("relTol", 0);
	picSolDict_.add("tolerance", 1.0e-6);
	picSolDict_.add("solver", "smoothSolver");
	picSolDict_.add("smoother", "symGaussSeidel");
}


bool pFlow::coupling::diffusion::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
	);

	auto& solidVol = solidVoldTmp.ref();
	
	size_t numPar = centerMass_.size();

	#pragma omp parallel for
	for(size_t i=0; i<numPar; i++)
	{
		const auto cellId = parCellIndex_[i];
		if( cellId >= 0 )
		{
			#pragma omp atomic
			solidVol[cellId] += 
				static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));
				
		}
	}

	auto picAlphaTmp = Foam::volScalarField::New(
		"picAlpha",
		this->mesh(),
		Foam::dimensioned("picAlpha", Foam::dimless, Foam::scalar(0)),
		"zeroGradient"
	);

	volScalarField& picAlpha = picAlphaTmp.ref();
	
	
	picAlpha.ref() = Foam::max( 1 - solidVol/this->mesh().V(), 0.0);
	picAlpha.correctBoundaryConditions();
	
	
	// start of Time loop
	for(Foam::label i=0; i<nSteps_; i++)
	{
		picAlpha.storeOldTime();
		fvScalarMatrix alphaEq
		(
			fvmDdt(picAlpha) - fvm::laplacian(DT_,picAlpha)
		);
		alphaEq.solve(picSolDict_);
	}

	this->ref() = picAlpha.internalField();

	return true;
}