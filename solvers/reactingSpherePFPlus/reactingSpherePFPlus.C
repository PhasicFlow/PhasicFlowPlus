/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
License
    This file is partially part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This file is partially part of phasicFlow.

    It is a free software for simulating granular and multiphase flows. You
    can redistribute it and/or modify it under the terms of GNU General
    Public License v3 or any other later versions.

    phasicFlow is distributed to help others in their research in the field of
    granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingSpherePFPlus

Description
    Transient incompressible solver for unresolved two-way coupling of
    spherical fluid-particle systems with heat exchange and multi-species
    heterogeneous reactions (chemical looping combustion). Extends
    unresolvedHeatSpherePFPlus with the transport of the gas species mass
    fractions, mixture properties computed from the local composition, and
    the multi-species reaction sources of the DEM particles.

    On the first PIMPLE outer iteration of every time step, the DEM side is
    advanced (contact, thermal and chemical calculations) and the heat and
    mass coupling sources are gathered, including the total gas mass source
    Sm and the molar production rate Sn of the gas. These sources are frozen
    for the remaining PIMPLE iterations, in which the momentum, pressure,
    mixture-property, species, and energy equations are solved in that
    order.

\*---------------------------------------------------------------------------*/

// OpenFOAM
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "CorrectPhi.H"

// phasicFlow
#include "unresolvedCouplingSystem.hpp"
#include "alphaTurbulentTransportModel.hpp"
#include "thermalMomentumSphereUnresolvedCouplingSystem.hpp"
#include "multiReactiveThermalMomentumSphereUnresolvedCouplingSystem.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    // Must precede the construction of the DEM system
    pFlow::Plus::processor::initMPI(argc, argv);

    // Creates the fields, the transport models and the coupling system
    #include "createFields.H"

    // Mixture properties of the initial composition and temperature
    #include "mixtureProperties.H"

    turbulence->validate();
    #include "CourantNo.H"

    // Reaction-rate timestep limit, read from the CLC dictionary of
    // system/fvSolution
    const dictionary& clcDict =
        mesh.solutionDict().subOrEmptyDict("CLC");

    const scalar maxReactionNo =
        clcDict.lookupOrDefault<scalar>("maxReactionNo", scalar(0.15));

    Info << "\n[CLC Solver] Reaction-rate timestep limit:\n"
         << "    maxReactionNo = " << maxReactionNo << "\n"
         << endl;

    Info << "\nStarting reactingSpherePFPlus Time Loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"

        #include "setDeltaT.H"

        // Reaction-rate-aware timestep constraint, applied to the time step
        // chosen by setDeltaT.H
        #include "ReactionNo.H"

        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        coupling.cfdTimers().start();

        // --- Pressure-velocity-species-temperature PIMPLE loop
        while (pimple.loop())
        {
            // DEM and coupling calculations are done only once per time
            // step, in the first PIMPLE iteration. The sources are frozen
            // for the subsequent iterations.
            if (pimple.firstIter())
            {
                coupling.cfdTimers().pause();

                const scalar t  = runTime.time().value();
                const scalar dt = runTime.deltaT().value();

                // Sync particle state
                coupling.getDataFromDEM(t, dt);

                // Map particles to Eulerian grid -> alpha
                coupling.calculatePorosity();

                alphaRhoPhi =
                    fvc::interpolate(U * coupling.alpha()) & mesh.Sf();

                // Compute drag + lift forces
                coupling.calculateMomentumCoupling();

                // Send previous step's fluid state to DEM GPU
                coupling.sendDataToDEM(t, dt);

                // Run DEM physics (contact, thermal, chemical)
                const bool demOk =
                    coupling.iterate
                    (
                        t, runTime.writeTime(), runTime.timeName()
                    );

                if (!demOk)
                {
                    FatalErrorIn(args.executable())
                        << "CRITICAL: DEM solver aborted. Halting OpenFOAM."
                        << exit(FatalError);
                }

                // Compute heat and mass coupling sources
                coupling.calculateHeatCoupling();
                coupling.calculateMassCoupling();

                // Total gas mass source: sum of the explicit source fields
                // of all species
                Sm = dimensionedScalar("zero", dimDensity/dimTime, scalar(0));

                const Foam::PtrList<Foam::volScalarField>& speciesSuFields =
                    coupling.speciesSu();

                forAll(speciesSuFields, i)
                {
                    Sm += speciesSuFields[i];
                }

                // Molar production rate of the gas: the species sources
                // divided by the molar masses
                Sn = dimensionedScalar
                (
                    "zero", dimMoles/dimVolume/dimTime, scalar(0)
                );

                forAll(gasSpeciesNames, j)
                {
                    const word SuName("Su_" + gasSpeciesNames[j]);

                    forAll(speciesSuFields, i)
                    {
                        if (speciesSuFields[i].name() == SuName)
                        {
                            Sn +=
                                speciesSuFields[i]
                              / dimensionedScalar
                                (
                                    "Mw", dimMass/dimMoles, gasMw[j]
                                );
                        }
                    }
                }

                coupling.cfdTimers().start();
            } // end pimple.firstIter()


            // Dynamic mesh update
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // Volume source of the gas released by the reactions, from the
            // molar volume of the current mixture
            Sv = Sn * Vm;

            // Momentum predictor
            #include "UEqn.H"

            // Pressure corrector
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            // Mixture properties (rho, Cp, kappa) from the current species
            // mass fractions, used by YEqn and TEqn
            #include "mixtureProperties.H"

            // Gas species mass fractions
            #include "YEqn.H"

            // Fluid energy equation
            #include "TEqn.H"

            // Turbulence model update
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }

        } // end pimple.loop()

        runTime.write();

        coupling.cfdTimers().end();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
             << "  ClockTime = "   << runTime.elapsedClockTime() << " s"
             << nl << endl;

    } // end runTime.run()

    Info << "End of reactingSpherePFPlus.\n" << endl;

    pFlow::Plus::processor::finalizeMPI();

    return 0;
}
