volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA * coupling.alpha())
  + MRF.zeroFilter(rAUf * fvc::interpolate(coupling.alpha())*fvc::ddtCorr(U, phi))
);

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU*coupling.alpha(), MRF);

volScalarField rAUAlpha2("rAUAlpha2",rAU * coupling.alpha() * coupling.alpha());

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUAlpha2, p) == fvc::div(phiHbyA) + fvc::ddt(coupling.alpha())
    );

    pEqn.setReference
    (
        pressureReference.refCell(),
        pressureReference.refValue()
    );

    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

U = fvc::reconstruct(phi)/coupling.alpha();
U.correctBoundaryConditions();
fvConstraints.constrain(U);
