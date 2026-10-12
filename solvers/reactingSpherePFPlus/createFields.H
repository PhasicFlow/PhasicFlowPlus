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

#include "createRDeltaT.H"

Info << "Reading basic fields: p, U, T\n" << endl;

volScalarField p
(
    IOobject
    (
        "p", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating multi-reactive thermal-momentum coupling system\n" << endl;

auto couplingPtr = pFlow::coupling::unresolvedCouplingSystem::create
(
    "multiReactiveSphere",
    "multiReactiveThermalMomentum",
    mesh,
    argc,
    argv
);

auto& coupling = couplingPtr();

// Typed pointer for the heat-coupling accessors (heatSp, heatSu)
auto* tCoupling = dynamic_cast
<
    pFlow::coupling::thermalMomentumSphereUnresolvedCouplingSystem*
>(&coupling);

if (!tCoupling)
{
    FatalErrorInFunction
        << "Coupling system '" << coupling.type()
        << "' is not a thermal-momentum coupling system."
        << exit(FatalError);
}

Info << "Calculating flux fields\n" << endl;

surfaceScalarField alphaRhoPhi
(
    IOobject
    (
        "alphaRhoPhi", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U * coupling.alpha()) & mesh.Sf()
);

surfaceScalarField phi
(
    IOobject
    (
        "phi", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf()
);

label  pRefCell  = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, alphaRhoPhi);

Info << "Reading physical properties from transportProperties\n" << endl;

dimensionedScalar rhoValue("rho", dimDensity, laminarTransport);
dimensionedScalar nuValue("nu", dimViscosity, laminarTransport);
dimensionedScalar CpValue
(
    "Cp", dimEnergy/dimMass/dimTemperature, laminarTransport
);
dimensionedScalar kappaValue
(
    "kappa", dimPower/dimLength/dimTemperature, laminarTransport
);
dimensionedScalar PrtValue("Prt", dimless, laminarTransport);

// Density, updated in place by mixtureProperties.H
volScalarField rho
(
    IOobject
    (
        "rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    rhoValue
);

// Specific heat capacity
volScalarField Cp
(
    IOobject
    (
        "Cp", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    CpValue
);

// Thermal conductivity
volScalarField kappa
(
    IOobject
    (
        "kappa", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    kappaValue
);

// Turbulent Prandtl number
volScalarField Prt
(
    IOobject
    (
        "Prt", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    PrtValue
);

dimensionedScalar Pr("Pr", (rhoValue * nuValue * CpValue) / kappaValue);

Info << "Solver Initialisation Summary:\n"
     << "    -> Calculated fluid Pr  : " << Pr.value() << "\n"
     << "    -> rho, Cp, kappa, Prt are registered fields.\n"
     << endl;

autoPtr<incompressible::alphaTurbulenceModel> turbulence
(
    incompressible::alphaTurbulenceModel::New
    (
        coupling.alpha(), U, alphaRhoPhi, phi, laminarTransport
    )
);

#include "createFvOptions.H"

autoPtr<surfaceVectorField> Uf;
if (mesh.dynamic())
{
    Uf = Foam::autoPtr<Foam::surfaceVectorField>::New
    (
        IOobject
        (
            "Uf", runTime.timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U * coupling.alpha())
    );
}

Info << "Reading gas species mass fractions from transportProperties\n"
     << endl;

wordList gasSpeciesNames(laminarTransport.lookup("gasSpecies"));
PtrList<volScalarField> Y(gasSpeciesNames.size());

forAll(gasSpeciesNames, i)
{
    Y.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Y." + gasSpeciesNames[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

// Per-species properties, in the order of gasSpecies
const List<scalar> gasMw(laminarTransport.lookup("gasMw"));
const List<scalar> gasCp(laminarTransport.lookup("gasCp"));
const List<scalar> gasKappa(laminarTransport.lookup("gasKappa"));

if (   gasMw.size()    != Y.size()
    || gasCp.size()    != Y.size()
    || gasKappa.size() != Y.size())
{
    FatalErrorInFunction
        << "gasMw, gasCp and gasKappa of transportProperties must have one "
        << "entry per gas species (" << Y.size() << ")."
        << exit(FatalError);
}

// Species obtained as 1 - sum(all others) instead of a transport equation
const word inertSpeciesName(laminarTransport.lookup("inertSpecies"));

label inertIndex = -1;
forAll(gasSpeciesNames, i)
{
    if (gasSpeciesNames[i] == inertSpeciesName)
    {
        inertIndex = i;
        break;
    }
}

if (inertIndex < 0)
{
    FatalErrorInFunction
        << "transportProperties/inertSpecies = '" << inertSpeciesName
        << "' is not in gasSpecies: " << gasSpeciesNames
        << exit(FatalError);
}

const dimensionedScalar Sc("Sc", dimless, laminarTransport);
const dimensionedScalar Sct("Sct", dimless, laminarTransport);

const dimensionedScalar pRef_abs("pRef", dimPressure, laminarTransport);

// Number of time steps in which YEqn.H rescaled the species fractions
label YscaleActivations = 0;

// Total gas mass source, assembled in reactingSpherePFPlus.C from the
// species sources of the particles
volScalarField Sm
(
    IOobject
    (
        "Sm", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sm", dimDensity/dimTime, scalar(0.0))
);

// Molar production rate of the gas [kmol/(m3.s)], assembled in
// reactingSpherePFPlus.C from the species sources of the particles
volScalarField Sn
(
    IOobject
    (
        "Sn", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Sn", dimMoles/dimVolume/dimTime, scalar(0.0))
);

// Molar volume of the gas mixture [m3/kmol], updated in mixtureProperties.H
volScalarField Vm
(
    IOobject
    (
        "Vm", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Vm", dimVolume/dimMoles, scalar(0.0))
);

// Volume source of the gas released by the reactions [1/s]: Sn*Vm
volScalarField Sv
(
    IOobject
    (
        "Sv", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sv", dimless/dimTime, scalar(0.0))
);

Info << "createFields.H complete.\n" << endl;
