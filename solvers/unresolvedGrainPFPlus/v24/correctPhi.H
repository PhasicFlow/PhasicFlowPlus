// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

correctUphiBCs(U, phi);

CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phi, U);
