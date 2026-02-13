volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA * coupling.alpha())
  + rAUf * fvc::interpolate(coupling.alpha())*fvc::ddtCorr(U, alphaRhoPhi, Uf) // ----> attention by Uf
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
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate((rAtU() - rAU)*coupling.alpha())*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

volScalarField rAUAlpha("rAUAlpha",rAtU() * coupling.alpha());

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUAlpha);

volScalarField rAUAlpha2("rAUAlpha2",rAtU() * coupling.alpha() * coupling.alpha());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUAlpha2, p) == fvc::div(phiHbyA) + fvc::ddt(coupling.alpha())  
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

//U = fvc::reconstruct(phi)/coupling.alpha();
U = HbyA - coupling.alpha() * rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// updates phi 
phi = linearInterpolate(U) & mesh.Sf();

// Correct Uf if the mesh is moving
if (mesh.dynamic())
{
    Uf() = fvc::interpolate(U * coupling.alpha());
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf() += n*(alphaRhoPhi/mesh.magSf() - (n & Uf()));
}

// Make the fluxes relative to the mesh motion
if (alphaRhoPhi.mesh().moving())
{
    alphaRhoPhi -= fvc::meshPhi(U * coupling.alpha());
}
