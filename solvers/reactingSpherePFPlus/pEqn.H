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

// Inverse of the diagonal coefficient of the momentum equation
volScalarField rAU(1.0 / UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

// Velocity without the pressure gradient (H operator over the diagonal)
volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

// Flux of HbyA
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA * coupling.alpha())
  + rAUf * fvc::interpolate(coupling.alpha())
  * fvc::ddtCorr(U, alphaRhoPhi, Uf)
);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0 / max(1.0 / rAU - UEqn.H1(), 0.1 / rAU);

    phiHbyA +=
        fvc::interpolate((rAtU() - rAU) * coupling.alpha())
      * fvc::snGrad(p) * mesh.magSf();

    HbyA -= (rAU - rAtU()) * fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

volScalarField rAUAlpha("rAUAlpha", rAtU() * coupling.alpha());

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUAlpha);

volScalarField rAUAlpha2
(
    "rAUAlpha2",
    rAtU() * coupling.alpha() * coupling.alpha()
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUAlpha2, p)
     ==
        fvc::div(phiHbyA)
      + fvc::ddt(coupling.alpha())
        // Volume source of the gas released by the reactions
      - (Sm / rho)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        alphaRhoPhi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - coupling.alpha() * rAU * fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Update the flux of the velocity
phi = linearInterpolate(U) & mesh.Sf();

// Correct Uf if the mesh is moving
if (mesh.dynamic())
{
    Uf.ref() = fvc::interpolate(U * coupling.alpha());
    surfaceVectorField n(mesh.Sf() / mesh.magSf());
    Uf.ref() += n * (alphaRhoPhi / mesh.magSf() - (n & Uf.ref()));
}

// Make the fluxes relative to the mesh motion
if (alphaRhoPhi.mesh().moving())
{
    alphaRhoPhi -= fvc::meshPhi(U * coupling.alpha());
}
