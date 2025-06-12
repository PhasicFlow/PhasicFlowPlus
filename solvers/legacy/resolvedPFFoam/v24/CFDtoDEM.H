PtrList<Foam::vector> particlesForce(nParticles);
PtrList<Foam::vector> particlesMoment(nParticles);

for(auto i=0; i<nParticles; ++i)
{

		// STL of the particle i
		const triSurface& STLi = particleSTLs[i];
		// Points of the STL of the particle i
		const pointField& STLpoints = STLi.points();

		// Center of mass of the particle i
		const point& centerOfMass = particlesCoM[i];

		// Fluid density
		dimensionedScalar rho("rho", dimDensity, laminarTransport);

		// Effective viscous stress
		tmp<Foam::volSymmTensorField> devRhoReff = -rho*turbulence->nu()*dev(twoSymm(fvc::grad(U)));

		// Pressure/viscous force on the particle
		Foam::vector particleiForce(0.0, 0.0, 0.0);

		// Pressure/viscous moment on the particle
		Foam::vector particleiMoment(0.0, 0.0, 0.0);

		forAll(STLi, facei)
		{
			const pointField& triPoints = STLi[facei].points(STLpoints);

			// Center of each triangle on the STL
			point faceiCentroid = (triPoints[0] + triPoints[1] + triPoints[2]) / 3.0;

			// Area vector of each triangle on the STL
			Foam::vector faceiArea = ((triPoints[1] - triPoints[0])^(triPoints[2] - triPoints[0])) / 2.0;

			// Cell ID with the center of triangle inside it
			label celli = mesh.findCell(faceiCentroid);

			// position of center of triangle relative to center of mass
			point centerToCoM = faceiCentroid - centerOfMass;

			// If the triangle center is in the CFD doamin
			if(celli != -1)
			{
				// Pressure force
				Foam::vector pressureForce = -rho.value() * faceiArea * p[celli]; //pf

				// Viscos force
				Foam::vector viscousForce = -faceiArea & devRhoReff().internalField()[celli]; // vf

				// Pressure moment
				Foam::vector pressureMoment = centerToCoM ^ pressureForce;

				// Viscous moment
				Foam::vector viscousMoment = centerToCoM ^ viscousForce;

				particleiForce += pressureForce + viscousForce;

				particleiMoment += pressureMoment + viscousMoment;
			}

		}

		reduce( particleiForce, sumOp<Foam::vector>() );
		reduce( particleiMoment, sumOp<Foam::vector>() );

		particlesForce.set(i, new Foam::vector(particleiForce));
		particlesMoment.set(i, new Foam::vector(particleiMoment));

		if (Pstream::master() || !Pstream::parRun())
		{
			auto& fluidForce = coupling.fluidForce();
			auto& fluidTorque = coupling.fluidTorque();

			fluidForce[i].x() = particlesForce[i][0];
			fluidForce[i].y() = particlesForce[i][1];
			fluidForce[i].z() = particlesForce[i][2];

			fluidTorque[i].x() = particlesMoment[i][0];
			fluidTorque[i].y() = particlesMoment[i][1];
			fluidTorque[i].z() = particlesMoment[i][2];
		}
}

coupling.sendFluidForceToDEM();
coupling.sendFluidTorqueToDEM();
