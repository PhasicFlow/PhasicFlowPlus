#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
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

Info<< "Reading field U\n" << endl;
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

// by default this is created on the master processor 
Info<< "Creating sphere coupling system for unresolved solver \n" << endl;
auto couplingPtr = pFlow::coupling::unresolvedCouplingSystem::create("grain", mesh, argc, argv);
auto& coupling = couplingPtr();

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U * coupling.alpha()) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, laminarTransport)
);

#include "createFvOptions.H"

autoPtr<surfaceVectorField> Uf;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf = Foam::autoPtr<Foam::surfaceVectorField>::New
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U * coupling.alpha())
    );
}
