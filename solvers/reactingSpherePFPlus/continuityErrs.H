/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of phasicFlow code. It is a free software for simulating
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions.

  phasicFlow is distributed to help others in their research in the field of
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

// Continuity error of the gas phase including the reaction mass source:
//
//     contErr = d(alpha)/dt + div(alpha*U) - Sm/rho
//
// Sm [kg/(m3.s)] is the total gas mass source of the reactions.

{
    volScalarField contErr
    (
        fvc::ddt(coupling.alpha())
      + fvc::div(alphaRhoPhi)
      - (Sm / rho)
    );

    // Sum of |contErr| weighted by cell volumes (local error metric)
    scalar sumLocalContErr =
        runTime.deltaTValue()
        * mag(contErr)().weightedAverage(mesh.V()).value();

    // Volume-weighted signed error (global conservation check)
    scalar globalContErr =
        runTime.deltaTValue()
        * contErr.weightedAverage(mesh.V()).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = "     << globalContErr
         << ", cumulative = " << cumulativeContErr
         << endl;
}
