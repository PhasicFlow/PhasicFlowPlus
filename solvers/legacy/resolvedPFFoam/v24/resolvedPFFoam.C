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
-------------------------------------------------------------------------------
 _____ _____ ____
|   __|   __|    \    Engineering    |
|   __|   __|  |  |   Fluid          |
|_____|__|  |____/    Dynamics       |  www.fluidDynamics.at

-------------------------------------------------------------------------------
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

    It is a free software for simulating
    granular and multiphase flows. You can redistribute it and/or modify it under
    the terms of GNU General Public License v3 or any other later versions.

    phasicFlow is distributed to help others in their research in the field of
    granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    resolvedSpherePFFoam

Description
    Transient incompressible solver for resolved simulation of spherical
    fluid-particle systems.

Author
    Bahram Haddadi, fluidDynamics.at.  All rights reserved - based on pimpleFoam

\*---------------------------------------------------------------------------*/

// OpenFOAM
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

// phasicFlow
#include "resolvedCouplingSystem.hpp"

// coupling
#include "STL.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    pFlow::Plus::processor::initMPI(argc, argv);

    #include "createFields.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;


    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Performing DEM calculation
        coupling.cfdTimers().start();
        coupling.iterate(runTime.time().value(), runTime.writeTime(), runTime.timeName());
        coupling.getDataFromDEM(runTime.time().value(), runTime.deltaT().value());

        // Transfering data from DEM to CFD
        #include "DEMtoCFD.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {

                mesh.update();

                coupling.cMesh().update(runTime.time().value(), runTime.deltaT().value());
                coupling.updateMeshBoxes();

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


            #include "UEqn.H"

            

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"

            } 

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // Transfering data from CFD to DEM
        #include "CFDtoDEM.H"

        if (runTime.write())
        {
          // Writing sphere STLs in a single STL file
            writeSTLs(particleSTLs, runTime.timeName());
        }

        coupling.cfdTimers().end();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    pFlow::Plus::processor::finalizeMPI();

    return 0;
}


// ************************************************************************* //
