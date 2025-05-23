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
	couplingSystem(shapeTypeName, mesh, argc, argv, true),
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

		const auto& pIDs = this->particleID();
		const auto& pCentersOfMass = this->centerMass();
		const auto& pLVelocities = this->particleVelocity();
		const auto& pRVelocities = this->particleRVelocity();

		// Creating the surface pointer for searching the cells inside the STL
		uniquePtr<Foam::triSurfaceSearch> querySurfPtr_ = makeUnique<triSurfaceSearch>(particleSTLs[i]);

		
		const auto& cellCenters = mesh.C();
		Foam::boolList insideCells = querySurfPtr_->calcInside(cellCenters);

		forAll(cellCenters, celli)
		{
			if(insideCells[celli])
			{
				//*** Alpha calculation needs to be updated to be between 0 and 1 ***//
				alpha[celli] = 0.0;
				particleID[celli] = pIDs[i];
				const Foam::vector centerOfMass {pCentersOfMass[i].x(), pCentersOfMass[i].y(), pCentersOfMass[i].z()};
				Foam::vector cellCenter {cellCenters[celli].x(), cellCenters[celli].y(), cellCenters[celli].z()};
				Foam::vector linearVelocity{pLVelocities[i].x(), pLVelocities[i].y(), pLVelocities[i].z()};
				Foam::vector rotVelocity{pRVelocities[i].x(), pRVelocities[i].y(), pRVelocities[i].z()};
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
	
	const auto& pCenterMass = this->centerMass();
	//Foam::PtrList<Foam::point> pCentersOfMass(this->particleCoMs());

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
				Foam::point centerToCoM = {
					faceiCentroid.x() - pCenterMass[i].x(),
					faceiCentroid.y() - pCenterMass[i].y(),
					faceiCentroid.z() - pCenterMass[i].z()};

				// If the triangle center is in the CFD domain
				if(celli != -1)
				{
					// Pressure force
					Foam::vector pressureForce = -rho[celli] * faceiArea * p[celli]; //pf

					// Viscous force
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

