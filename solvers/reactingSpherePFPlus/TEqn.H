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

// Energy equation of the gas phase

{
    // Effective thermal diffusivity [m2/s]: kappa/(rho*Cp) + nut/Prt
    volScalarField alphaEff
    (
        "alphaEff",
        (kappa / (rho * Cp)) + turbulence->nut() / Prt
    );

    const volScalarField& alpha = coupling.alpha();

    fvScalarMatrix TEqn
    (
        fvm::ddt(alpha, T)
      + fvm::div(alphaRhoPhi, T)
      - fvm::laplacian(alpha * alphaEff, T)
    );

    // Particle-fluid convection. heatSp() is negative (-h*Ap/V), so the
    // diagonal coefficient -heatSp/(rho*Cp) is positive.
    TEqn += fvm::Sp(-(tCoupling->heatSp() / (rho * Cp)), T)
          - (tCoupling->heatSu() / (rho * Cp));

    // The gas released by the reactions enters at the local fluid
    // temperature, so the term T*Sv contained in the matrix above is removed
    TEqn -= Sv * T;

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);
}
