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

/**
 * @file createFields.H
 * @brief Initializes all primary CFD fields and the DEM coupling system
 * for the unresolved thermal (non-reactive) CFD-DEM solver.
 */

#include "createRDeltaT.H"

// =========================================================================//
// Section 1: Primary Field Instantiation (p, U, T)
// =========================================================================//
Info << "Reading basic fields: p, U, T\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// =========================================================================//
// Section 2: Coupling System Instantiation
// =========================================================================//
Info << "Creating thermal-momentum coupling system\n" << endl;

auto couplingPtr = pFlow::coupling::unresolvedCouplingSystem::create
(
    "thermalSphere",
    "thermalMomentum",
    mesh,
    argc,
    argv
);

auto& coupling = couplingPtr();

// Extract the typed pointer for heat-coupling accessors.
// Forces a fatal error if the wrong coupling type is configured to
// prevent the solver from silently failing during the energy update.
using ThermalCoupling =
    pFlow::coupling::thermalMomentumSphereUnresolvedCouplingSystem;

auto* tCoupling = dynamic_cast<ThermalCoupling*>(&coupling);

if (!tCoupling)
{
    FatalErrorInFunction
        << "Coupling system '" << coupling.type()
        << "' does not implement the thermal coupling interface.\n"
        << "Ensure shape type = 'thermalSphere' and coupling type "
        << "= 'thermalMomentum'."
        << exit(FatalError);
}

// =========================================================================//
// Section 3: Volumetric Flux Definitions
// =========================================================================//
Info << "Reading/calculating face volumetric flux field alphaRhoPhi\n"
     << endl;

// rho is considered to be one for the flux definition below
surfaceScalarField alphaRhoPhi
(
    IOobject
    (
        "alphaRhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U * coupling.alpha()) & mesh.Sf()
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

label  pRefCell  = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

// =========================================================================//
// Section 4: Transport & Thermophysical Properties
// =========================================================================//
singlePhaseTransportModel laminarTransport(U, alphaRhoPhi);

Info << "Reading/Registering physical properties in ObjectRegistry\n"
     << endl;

dimensionedScalar rhoValue
(
    "rho", dimDensity, laminarTransport
);

dimensionedScalar nuValue
(
    "nu", dimViscosity, laminarTransport
);

dimensionedScalar CpValue
(
    "Cp", dimEnergy/dimMass/dimTemperature, laminarTransport
);

dimensionedScalar kappaValue
(
    "kappa", dimPower/dimLength/dimTemperature, laminarTransport
);

dimensionedScalar PrtValue
(
    "Prt", dimless, laminarTransport
);

// Density field
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

// Kinematic viscosity
volScalarField nu
(
    IOobject
    (
        "nu", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    nuValue
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
     << "    -> Calculated fluid Pr  : " << Pr.value()    << "\n"
     << "    -> nu, Cp, kappa, Prt are registered fields.\n"
     << endl;

// =========================================================================//
// Section 5: Turbulence Model
// =========================================================================//
autoPtr<incompressible::alphaTurbulenceModel> turbulence
(
    incompressible::alphaTurbulenceModel::New
    (
        coupling.alpha(), U, alphaRhoPhi, phi, laminarTransport
    )
);

#include "createFvOptions.H"

// =========================================================================//
// Section 6: Face-Velocity Field (Dynamic Mesh Only)
// =========================================================================//
autoPtr<surfaceVectorField> Uf;
if (mesh.dynamic())
{
    Info << "Constructing face velocity Uf\n" << endl;

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

Info << "createFields.H complete.\n" << endl;


