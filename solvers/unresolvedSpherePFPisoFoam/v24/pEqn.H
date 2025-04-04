volScalarField rAU(1.0/UEqn.A());
volScalarField rAUAlpha("rAUAlpha",rAU * coupling.alpha());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA * coupling.alpha())
  + rAUf * fvc::interpolate(coupling.alpha())*fvc::ddtCorr(U, phi)
);

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUAlpha);

volScalarField rAUAlpha2("rAUAlpha2",rAU * coupling.alpha() * coupling.alpha());

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUAlpha2, p) == fvc::div(phiHbyA) + fvc::ddt(coupling.alpha())
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

U = fvc::reconstruct(phi)/coupling.alpha();
U.correctBoundaryConditions();
fvOptions.correct(U);
