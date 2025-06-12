// Reading data from DEM to CFD on the master processor
int nParticles = 0;
int resolution = 20; // this should also be defined in a dictionary as a relative parameter

volScalarField& alpha = coupling.alpha();

if (Pstream::master() || !Pstream::parRun())
{
	nParticles = coupling.numParticles();
}

typedef DynamicList<scalar> listScalar;
typedef DynamicList<point> listPoint;
typedef DynamicList<triSurface> listTriSurface;

List<listScalar> pDiameters(nParticles);
List<listPoint> pCenterOfMasses(nParticles);
List<listPoint> pLinVelocities(nParticles);
//Needs to be imported from DEM side
List<listPoint> pRotVelocities(nParticles);
List<listTriSurface> sphereSTLs(nParticles);

auto particleDiameters = coupling.particleDiameter();
auto particleCenterOfMassess = coupling.centerMass();
auto particleVelocities = coupling.particleVelocity();

if (Pstream::master() || !Pstream::parRun())
{
	for(auto i=0; i<nParticles; ++i)
	{
      // Transfers the particle diameters to OpenFOAM
      pDiameters[i].append(particleDiameters[i]);

      // Transfers the particle centers of mass to OpenFOAM
      point centerOfMass;
      {
				Foam::vector& point(centerOfMass);
				point[0] = particleCenterOfMassess[i].x();
				point[1] = particleCenterOfMassess[i].y();
				point[2] = particleCenterOfMassess[i].z();
			}
      pCenterOfMasses[i].append(centerOfMass);

      // Transfers the particle linear velocities to OpenFOAM
      point lVelocity;
      {
				Foam::vector& point(lVelocity);
				point[0] = particleVelocities[i].x();
				point[1] = particleVelocities[i].y();
				point[2] = particleVelocities[i].z();
			}
      pLinVelocities[i].append(lVelocity);

      // Transfers the particle rotational velocities to OpenFOAM
      point rVelocity;
      {
				Foam::vector& point(rVelocity);
				point[0] = 0.0; // Dummy value for checking the code
				point[1] = 0.0;
				point[2] = 0.0;
			}
      pRotVelocities[i].append(rVelocity);

      // Generating sphere STLs
      triSurface sphereSTL
    = generateSphere
      (
        pCenterOfMasses[i][0],
        pDiameters[i][0],
        resolution
      );
      sphereSTLs[i].append(sphereSTL);
	}
}

List<List<listPoint > > shareAllCenterOfMass(Pstream::nProcs()) ;
List<List<listPoint > > shareAllLVelocities(Pstream::nProcs()) ;
List<List<listPoint > > shareAllRVelocities(Pstream::nProcs()) ;
List<List<listTriSurface> >  shareAllSTLs(Pstream::nProcs()) ;

shareAllCenterOfMass[Pstream::myProcNo()].transfer(pCenterOfMasses);
shareAllLVelocities[Pstream::myProcNo()].transfer(pLinVelocities);
shareAllRVelocities[Pstream::myProcNo()].transfer(pRotVelocities);
shareAllSTLs[Pstream::myProcNo()].transfer(sphereSTLs);

Pstream::gatherList(shareAllCenterOfMass);
Pstream::gatherList(shareAllLVelocities);
Pstream::gatherList(shareAllRVelocities);
Pstream::gatherList(shareAllSTLs);

Pstream::scatterList(shareAllCenterOfMass);
Pstream::scatterList(shareAllLVelocities);
Pstream::scatterList(shareAllRVelocities);
Pstream::scatterList(shareAllSTLs);

PtrList<point> particlesCoM(nParticles);
PtrList<point> particlesLVel(nParticles);
PtrList<point> particlesRVel(nParticles);
PtrList<triSurface> particleSTLs(nParticles);

for(auto i=0; i<nParticles; ++i)
{
    	particlesCoM.set(i, new point(shareAllCenterOfMass[0][i][0]));
    	particlesLVel.set(i, new point(shareAllLVelocities[0][i][0]));
    	particlesRVel.set(i, new point(shareAllRVelocities[0][i][0]));
      particleSTLs.set(i, new triSurface(shareAllSTLs[0][i][0]));
}

// Initializing alpha and particleID fields
forAll(alpha, celli)
{
	alpha[celli] = 1.0;
	particleID[celli] = 0.0;
	Us[celli] = 0.0*Us[celli];
}

// Updating alpha, particleID and solid velocity field based on the DEM particles
for(auto i=0; i<nParticles; ++i)
{
	const triSurfaceSearch* querySurfPtr_(new triSurfaceSearch(particleSTLs[i]));

	boolList insideCells = querySurfPtr_->calcInside(mesh.C());

	forAll(mesh.C(), celli)
	{
		if(insideCells[celli])
		{
			alpha[celli] = 0.0;
			particleID[celli] = i+1;
			Foam::vector centerOfMass = particlesCoM[i];
			Foam::vector cellCenter = mesh.C()[celli];
			Foam::vector linearVelocity = particlesLVel[i];
			Foam::vector rotVelocity = particlesRVel[i];
			Us[celli] = linearVelocity + (rotVelocity^(cellCenter-centerOfMass));
		}
	}
}

alpha.correctBoundaryConditions();
particleID.correctBoundaryConditions();
Us.correctBoundaryConditions();
