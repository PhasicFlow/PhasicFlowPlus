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
#include "kinematicMomentumTransportModel.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"

// phasicFlow
#include "couplingSystem.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        coupling.iterate(runTime.time().value(), 1.0, runTime.timeName());
        coupling.getDaraFromDEM();
        coupling.calculatePorosity();
        volScalarField& voidFraction = coupling.alpha();


        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            fvModels.correct();

            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();
        coupling.collectFluidForce();
        coupling.collectFluidTorque();
        // For GPU parallelization
        coupling.sendFluidForceToDEM();
        coupling.sendFluidTorqueToDEM();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
