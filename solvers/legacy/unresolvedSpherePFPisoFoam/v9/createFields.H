
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        Foam::timeName(runTime),
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
        Foam::timeName(runTime),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// by default this is created on the master processor 
Info<< "Creating sphere coupling system for unresolved solver \n" << endl;
auto couplingPtr = pFlow::coupling::unresolvedCouplingSystem::create("sphere", mesh, argc, argv);
auto& coupling = couplingPtr();

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        Foam::timeName(runTime),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U * coupling.alpha()) & mesh.Sf()
);


pressureReference pressureReference(p, piso.dict());

mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        Foam::timeName(runTime),
        mesh
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, laminarTransport.lookup("rho"))
);

#include "createFvModels.H"
#include "createFvConstraints.H"


