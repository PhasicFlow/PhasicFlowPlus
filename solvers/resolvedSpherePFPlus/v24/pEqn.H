volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    volScalarField ddtAlpha = fvc::ddt(alpha);
    // Calculation of new ddtAlpha by using divergence of particle fluxes instead
    // volScalarField ddtAlpha=fvc::div(Us*(1.0-alpha));

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) + ddtAlpha
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(p.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();

// Correction for DEM IBM
U=(1.0-alpha)*Us+alpha*U;

fvScalarMatrix IBVolFractionEqn
(
    fvm::laplacian(IBVolFraction)
 ==
    fvc::ddt(alpha)
  + fvc::div(U)
);

if(IBVolFraction.needReference())
{
    setRefCell(IBVolFraction, pimple.dict(), pRefCell, pRefValue);
}

IBVolFractionEqn.solve();

U=U-fvc::grad(IBVolFraction);
U.correctBoundaryConditions();

p=p+IBVolFraction/U.mesh().time().deltaT();
p.correctBoundaryConditions();

fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

