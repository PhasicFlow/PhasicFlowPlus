// Solid->Fluid interaction: Just the particles that are inside each CFD sub-domain are considered

volScalarField& alpha = coupling.alpha();

PtrList<scalar> pDiameters(coupling.particleDiameters());
PtrList<point> pCentersOfMass(coupling.particleCoMs());

// Creating sphere STLs for each particle
PtrList<triSurface> particleSTLs(coupling.numParticles());

forAll(pDiameters, i)
{
	triSurface sphereSTL 
  = 
    generateSphere
	(
		pCentersOfMass[i],
		pDiameters[i],
		readScalar(coupling.lookup("sphereSTLResolution"))
	);
	particleSTLs.set(i, new triSurface(sphereSTL));
}

// Calculating the solid to fluid interaction
coupling.calculateSolidInteraction
(
	mesh,
	particleSTLs,
	alpha,
	particleID,
	Us
);

