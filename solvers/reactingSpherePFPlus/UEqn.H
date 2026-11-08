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

// Momentum equation of the gas phase

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(alphaRhoPhi, U)
  + turbulence->divDevReff(U)

    // Implicit drag: Sp >= 0 adds Sp/rho to the matrix diagonal
  + fvm::Sp(coupling.Sp() / rho, U)

 ==

    // Explicit drag source (includes lift force components)
  - coupling.Su() / rho

    // Gas released by the reactions carries the local fluid velocity
  + Sv * U

  + fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -coupling.alpha() * fvc::grad(p));

    fvOptions.correct(U);
}
