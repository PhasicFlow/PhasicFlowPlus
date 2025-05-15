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

Info<< "Creating/Reading field Us\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Us", dimVelocity, vector(0.0,0.0,0.0))
);

Info<< "Reading field immersed boundary volume fraction\n" << endl;
volScalarField IBVolFraction
(
    IOobject
    (
        "IBVolFraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating/Reading field particleID\n" << endl;
volScalarField particleID
(
    IOobject
    (
        "particleID",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("particleID", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

Info<< "Creating sphere coupling system for resolved solver \n" << endl;
pFlow::coupling::resolvedCouplingSystem coupling("sphere", mesh, argc, argv);

# include "createPhi.H"

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
    dimensionedScalar
    (
        "rho", 
        dimDensity, 
        readScalar(laminarTransport.lookup("rho"))
    )
);

#include "createFvOptions.H"

autoPtr<surfaceVectorField> Uf;

if (mesh.dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf.reset
    (
        new surfaceVectorField
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
        ) 
    );
}
