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

    It is a free software for simulating 
    granular and multiphase flows. You can redistribute it and/or modify it under
    the terms of GNU General Public License v3 or any other later versions. 
 
    phasicFlow is distributed to help others in their research in the field of 
    granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    unresolvedSpherePFFoam

Description
    Transient incompressible solver for unresolved simulation of spherical 
    fluid-particle systems.


\*---------------------------------------------------------------------------*/

// OpenFOAM
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// phasicFlow
#include "unresolvedCouplingSystem.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    pFlow::Plus::processor::initMPI(argc, argv);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        auto t = runTime.time().value();
        auto dt = runTime.deltaT().value();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //coupling.cMesh().update(t, dt);
        
        coupling.getDataFromDEM(t, dt);
        coupling.calculatePorosity();
        coupling.calculateFluidInteraction();
        coupling.sendDataToDEM(t, dt);

        Info<<"Iterating DEM up to time "<< runTime.time().value()<<endl;
        coupling.iterate(t, runTime.writeTime(), runTime.timeName());

        coupling.cfdTimers().start();
        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {

            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

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
