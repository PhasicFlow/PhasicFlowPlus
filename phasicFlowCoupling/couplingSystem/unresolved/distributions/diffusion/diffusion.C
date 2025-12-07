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
#include "couplingMesh.hpp"


Foam::tmp<Foam::fvMatrix<Foam::scalar>> pFlow::coupling::diffusion::fvmDdt
(
    const Foam::volScalarField& sField
)const
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

Foam::tmp<Foam::fvMatrix<Foam::vector>> pFlow::coupling::diffusion::fvmDdt
(
    const Foam::volVectorField& sField
)const
{

    Foam::tmp<Foam::fvMatrix<Foam::vector>> tfvm
    (
        new Foam::fvMatrix<Foam::vector>
        (
            sField,
            sField.dimensions()*Foam::dimVol/Foam::dimTime
        )
    );

    Foam::fvMatrix<Foam::vector>& fvm = tfvm.ref();

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

void pFlow::coupling::diffusion::constructLists
(
    const Foam::scalar searchLen, 
    const Foam::label maxLayers
)
{
}

void pFlow::coupling::diffusion::constructLists(const Foam::label maxLayers)
{
}

pFlow::coupling::diffusion::diffusion
(
    const Foam::dictionary &parrentDict,
    const couplingMesh &cMesh,
    const Plus::centerMassField &centerMass
)
: 
    distributionBase
    (
        false, 
        parrentDict, 
        cMesh, 
        centerMass
    ),
    nSteps_
    (
        Foam::max(1, parrentDict.subDict("diffusionInfo").template get<Foam::label>("nSteps"))
    ),
    standardDeviation_
    (
        parrentDict.subDict("diffusionInfo").template get<Foam::scalar>("standardDeviation")
    ),
    intTime_(1.0),
    dt_
    (
        "dt", 
        Foam::dimTime, 
        intTime_ / nSteps_
    ),
    DT_
    (
        "DT", 
        Foam::dimDynamicViscosity / Foam::dimDensity, 
        standardDeviation_ * standardDeviation_ / 4
    ),
    smoothSolDict_("smoothSolDict")
{
    smoothSolDict_.add("relTol", 0);
	smoothSolDict_.add("tolerance", 1.0e-7);
	smoothSolDict_.add("solver", "smoothSolver");
	smoothSolDict_.add("smoother", "symGaussSeidel");
    int log = parrentDict.subDict("diffusionInfo").lookupOrDefault<int>("log", 0);
    smoothSolDict_.add("log", log);
}

void pFlow::coupling::diffusion::updateWeights(const Plus::procCMField<real> &parDiameter)
{
}

void pFlow::coupling::diffusion::smoothenField(Foam::volScalarField &field)const
{
    auto fldName = Foam::IOobject::groupName(field.name(), "smooth");
    auto tmpSmooth = Foam::volScalarField::New(
		fldName,
		field.mesh(),
		Foam::dimensioned(fldName, field.dimensions(), Foam::scalar(0)),
		"zeroGradient"
	);

	Foam::volScalarField& smooth = tmpSmooth.ref();
	
    forAll(smooth, celli)
    {
        smooth[celli] =  field[celli];
    }
	
    smooth.correctBoundaryConditions();
	
	
	// start of Time loop
    Foam::Info<< Blue_Text("Diffusion: smoothing field ") << 
                 Blue_Text(field.name()) << 
                 Blue_Text(" with ") << 
                 Yellow_Text(nSteps_ )<< 
                 Blue_Text(" steps.") << Foam::endl;
                
	for(Foam::label i=0; i<nSteps_; i++)
	{
		smooth.storeOldTime();
		Foam::fvScalarMatrix smoothEq
		(
			fvmDdt(smooth) - Foam::fvm::laplacian(DT_,smooth)
		);
		smoothEq.solve(smoothSolDict_);
	}

    forAll(smooth, celli)
    {
        field[celli] = smooth[celli];
    }
}

void pFlow::coupling::diffusion::smoothenField(Foam::volVectorField& field)const
{
    auto fldName = Foam::IOobject::groupName(field.name(), "smooth");
    auto tmpSmooth = Foam::volVectorField::New(
		fldName,
		field.mesh(),
		Foam::dimensioned(fldName, field.dimensions(), Foam::vector(0,0,0)),
		"zeroGradient"
	);

	Foam::volVectorField& smooth = tmpSmooth.ref();
	
    forAll(smooth, celli)
    {
        smooth[celli] =  field[celli];
    }
		
    smooth.correctBoundaryConditions();
	
	
	// start of Time loop
    Foam::Info<< Blue_Text("Diffusion: smoothing field ") << 
                Blue_Text(field.name()) << 
                Blue_Text(" with ") << 
                Yellow_Text(nSteps_ )<< 
                Blue_Text(" steps.") << Foam::endl;

	for(Foam::label i=0; i<nSteps_; i++)
	{
		smooth.storeOldTime();
		Foam::fvVectorMatrix smoothEq
		(
			fvmDdt(smooth) - Foam::fvm::laplacian(DT_,smooth)
		);
		smoothEq.solve(smoothSolDict_);
	}

    forAll(smooth, celli)
    {
        field[celli] = smooth[celli];
    }

}
