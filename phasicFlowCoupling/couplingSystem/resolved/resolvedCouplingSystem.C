/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
 _____ _____ ____
|   __|   __|    \    Engineering    |
|   __|   __|  |  |   Fluid          |
|_____|__|  |____/    Dynamics       |  www.fluidDynamics.at

-------------------------------------------------------------------------------
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#include "resolvedCouplingSystem.hpp"

pFlow::coupling::resolvedCouplingSystem::resolvedCouplingSystem
(
	word shapeTypeName, 
	Foam::fvMesh& mesh,
	int argc, 
	char* argv[]
)
:
	couplingSystem(shapeTypeName, mesh, argc, argv),
	alpha_
	(
		Foam::IOobject
	    (
	        "alpha",
	        mesh.time().timeName(),
	        mesh,
	        Foam::IOobject::MUST_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	mesh
	)
{}


Foam::PtrList<Foam::scalar> 
pFlow::coupling::resolvedCouplingSystem::particleIDs()
{
	int nParticles = this->numParticles();
	Foam::PtrList<Foam::scalar> pIDs(nParticles);

	for(auto i=0; i<nParticles; ++i)
	{
		pIDs.set(i, new scalar(this->couplingSystem::particleID()[i]));
	}

	return pIDs;
}

Foam::PtrList<Foam::scalar> 
pFlow::coupling::resolvedCouplingSystem::particleDiameters()
{
	int nParticles = this->numParticles();
	Foam::PtrList<Foam::scalar> pDiameters(nParticles);

	for(auto i=0; i<nParticles; ++i)
	{
		pDiameters.set(i, new scalar(this->couplingSystem::particleDiameter()[i]));
	}

	return pDiameters;
}

Foam::PtrList<Foam::point> 
pFlow::coupling::resolvedCouplingSystem::particleCoMs()
{
	int nParticles = this->numParticles();
	Foam::PtrList<Foam::point> pCentersOfMass(nParticles);

	for(auto i=0; i<nParticles; ++i)
	{
		Foam::point centerOfMass;
		{
			Foam::vector& point(centerOfMass);
			point[0] = this->couplingSystem::centerMass()[i].x();
			point[1] = this->couplingSystem::centerMass()[i].y();
			point[2] = this->couplingSystem::centerMass()[i].z();
		}
		pCentersOfMass.set(i, new point(centerOfMass));
	}

	return pCentersOfMass;
}

Foam::PtrList<Foam::point> 
pFlow::coupling::resolvedCouplingSystem::particleLVelocities()
{
	int nParticles = this->numParticles();
	Foam::PtrList<Foam::point> pLVelocities(nParticles);

	for(auto i=0; i<nParticles; ++i)
	{
		Foam::point lVelocity;
		{
			Foam::vector& point(lVelocity);
			point[0] = this->couplingSystem::particleVelocity()[i].x();
			point[1] = this->couplingSystem::particleVelocity()[i].y();
			point[2] = this->couplingSystem::particleVelocity()[i].z();
		}
		pLVelocities.set(i, new point(lVelocity));
	}

	return pLVelocities;
}

Foam::PtrList<Foam::point>
pFlow::coupling::resolvedCouplingSystem::particleRVelocities()
{
	int nParticles = this->numParticles();
	Foam::PtrList<Foam::point> pRVelocities(nParticles);

	for(auto i=0; i<nParticles; ++i)
	{
		Foam::point rVelocity;
		{
			Foam::vector& point(rVelocity);
			point[0] = this->couplingSystem::particleRVelocity()[i].x();
			point[1] = this->couplingSystem::particleRVelocity()[i].y();
			point[2] = this->couplingSystem::particleRVelocity()[i].z();
		}
		pRVelocities.set(i, new point(rVelocity));
	}

	return pRVelocities;
}

void pFlow::coupling::resolvedCouplingSystem::calculateSolidInteraction
(
	Foam::fvMesh& mesh,
	Foam::PtrList<Foam::triSurface>& particleSTLs,
	Foam::volScalarField& alpha,
	Foam::volScalarField& particleID,
	Foam::volVectorField& Us
) 
{
	int nParticles = this->numParticles();

	// Initializing alpha and particleID fields
	forAll(alpha, celli)
	{
		alpha[celli] = 1.0;
		particleID[celli] = 0.0;
		Us[celli] = 0.0*Us[celli];
	}

	// Updating CFD fields for each particle
	for(auto i=0; i<nParticles; ++i)
	{
		// Converting PhasicFlow fields to OpenFOAM fields
		Foam::PtrList<Foam::scalar> pIDs(this->particleIDs());
		Foam::PtrList<Foam::point> pCentersOfMass(this->particleCoMs());
		Foam::PtrList<Foam::point> pLVelocities(this->particleLVelocities());
		Foam::PtrList<Foam::point> pRVelocities(this->particleRVelocities());

		// Creating the surface pointer for searching the cells inside the STL
		const triSurfaceSearch* querySurfPtr_(new triSurfaceSearch(particleSTLs[i]));

		Foam::boolList insideCells = querySurfPtr_->calcInside(mesh.C());

		forAll(mesh.C(), celli)
		{
			if(insideCells[celli])
			{
				//*** Alpha calculation needs to be updated to be between 0 and 1 ***//
				alpha[celli] = 0.0;
				particleID[celli] = pIDs[i];
				Foam::vector centerOfMass = pCentersOfMass[i];
				Foam::vector cellCenter = mesh.C()[celli];
				Foam::vector linearVelocity = pLVelocities[i];
				Foam::vector rotVelocity = pRVelocities[i];
				Us[celli] = linearVelocity + (rotVelocity^(cellCenter-centerOfMass));
			}
		}
	}

	// Correcting and updating the boundary conditions
	alpha.correctBoundaryConditions();
	particleID.correctBoundaryConditions();
	Us.correctBoundaryConditions();
}

void pFlow::coupling::resolvedCouplingSystem::calculateFluidInteraction
(
	Foam::volScalarField& p,
	Foam::volScalarField& rho,
	Foam::volSymmTensorField& devRhoReff,
	Foam::PtrList<Foam::triSurface>& particleSTLs
) 
{
	int nParticles = this->numParticles();

	Foam::PtrList<Foam::point> pCentersOfMass(this->particleCoMs());

	for(auto i=0; i<nParticles; ++i)
	{
			// STL of the particle i
			const Foam::triSurface& STLi = particleSTLs[i];

			// Pressure/viscous force on the particle
			Foam::vector particleiForce(0.0, 0.0, 0.0);

			// Pressure/viscous moment on the particle
			Foam::vector particleiMoment(0.0, 0.0, 0.0);

			Foam::label celli = 0;

			forAll(STLi, facei)
			{
				const Foam::pointField& triPoints = STLi[facei].points(STLi.points());

				// Center of each triangle on the STL
				Foam::point faceiCentroid = (triPoints[0] + triPoints[1] + triPoints[2]) / 3.0;

				// Area vector of each triangle on the STL
				Foam::vector faceiArea = ((triPoints[1] - triPoints[0])^(triPoints[2] - triPoints[0])) / 2.0;

				// Cell ID with the center of triangle inside it
				celli = this->cMesh().findPointInCellTree(faceiCentroid, celli);

				// position of center of triangle relative to center of mass
				Foam::point centerToCoM = faceiCentroid - pCentersOfMass[i];

				// If the triangle center is in the CFD doamin
				if(celli != -1)
				{
					// Pressure force
					Foam::vector pressureForce = -rho[celli] * faceiArea * p[celli]; //pf

					// Viscos force
					Foam::vector viscousForce = -faceiArea & devRhoReff[celli]; // vf

					// Pressure moment
					Foam::vector pressureMoment = centerToCoM ^ pressureForce;

					// Viscous moment
					Foam::vector viscousMoment = centerToCoM ^ viscousForce;

					particleiForce += pressureForce + viscousForce;

					particleiMoment += pressureMoment + viscousMoment;
				}

			}

			// Sending fluid force on particle i to PhsicFlow
			this->fluidForce()[i].x() = particleiForce[0];
			this->fluidForce()[i].y() = particleiForce[1];
			this->fluidForce()[i].z() = particleiForce[2];

			// Sending fluid torque on particle i to PhsicFlow
			this->fluidTorque()[i].x() = particleiMoment[0];
			this->fluidTorque()[i].y() = particleiMoment[1];
			this->fluidTorque()[i].z() = particleiMoment[2];
	}
}

