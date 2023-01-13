
// by default this is created on the master processor 
Info<< "Creating sphere coupling system for unresolved solver \n" << endl;
pFlow::coupling::couplingSystem coupling("sphereDEMSystem", mesh, argc, argv);


//#include "createRDeltaT.H"

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

//#include "createPhi.H"

pressureReference pressureReference(p, pimple.dict());

//mesh.setFluxRequired(p.name());


//singlePhaseTransportModel laminarTransport(U, phi);

/*autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);*/


//#include "createFvModels.H"
//#include "createFvConstraints.H"
