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

/**
 * @file UEqn.H
 * @brief Solves the fluid-phase momentum equation for the unresolved
 * thermal (non-reactive) CFD-DEM solver.
 */

// Solve the Momentum equation
volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(alphaRhoPhi, U)
  //- fvc::div(alphaNuEff * dev(fvc::grad(U)().T()))
  + turbulence->divDevReff(U)
  + fvm::Sp(coupling.Sp()/rho, U)
 ==
  - coupling.Su()/rho
  + fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == - coupling.alpha() * fvc::grad(p));

    fvOptions.correct(U);
}


