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
 * @file TEqn.H
 * @brief Fluid-phase energy (temperature) equation for the unresolved
 * thermal (non-reactive) CFD-DEM solver.
 *
 * @details
 * Formulates and solves the macroscopic energy balance for the fluid
 * phase, integrating implicit/explicit convective coupling from the DEM
 * particles.
 *
 * Energy Source Term — Particle-fluid convection (Implicit + Explicit):
 * `fvm::Sp(-(heatSp)/(rho*Cp), T) - heatSu/(rho*Cp)`
 * where:
 * - heatSp = -h*Ap/V [W/(m³·K)] (Negative scalar)
 * - heatSu =  h*Ap*Tp/V [W/m³] (Positive scalar when Tp > 0)
 */

{
    //--- effective thermal diffusivity --------------------------------------
    // α_eff = κ/(ρ·Cp) + ν_t/Prt   [m²/s]
    // First term:  molecular (laminar) thermal diffusivity
    // Second term: turbulent thermal diffusivity
    volScalarField alphaEff
    (
        "alphaEff",
        (kappa / (rho * Cp)) + turbulence->nut() / Prt
    );

    const volScalarField& alpha = coupling.alpha();

    //--- base fluid-phase energy equation matrix ----------------------------
    // ∂(α·T)/∂t + ∇·(α·U·T) − ∇·(α·α_eff·∇T) = Sources
    fvScalarMatrix TEqn
    (
        fvm::ddt(alpha, T)
      + fvm::div(alphaRhoPhi, T)
      - fvm::laplacian(alpha * alphaEff, T)
    );

    //--- particle-fluid convective heat exchange ----------------------------
    // Numerical stability note: `heatSp()` is evaluated and stored in the
    // coupling layer as a NEGATIVE quantity (e.g., -hAp/V). OpenFOAM's
    // `fvm::Sp(coeff, T)` ADDS the coefficient to the matrix diagonal. To
    // ensure diagonal dominance and solver stability, the diagonal
    // coefficient must be strictly POSITIVE. Therefore, a minus sign is
    // applied: `fvm::Sp(-(heatSp)/..., T)`.
    TEqn += fvm::Sp(-(tCoupling->heatSp() / (rho * Cp)), T)
          - (tCoupling->heatSu() / (rho * Cp));

    //--- matrix relaxation, constraint, and solve ---------------------------
    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);
}
