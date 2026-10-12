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

// Mixture properties from the local composition and temperature: ideal-gas
// density, mass-fraction-weighted heat capacity and thermal conductivity.
// The kinematic viscosity keeps its transportProperties value.

{
    // Mixture molecular weight [kg/kmol]: 1/Mw_mix = sum(Y_i/Mw_i)
    volScalarField invMwMix
    (
        IOobject
        (
            "invMwMix", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMoles / dimMass, 0.0)
    );

    forAll(Y, i)
    {
        const dimensionedScalar Mw_i
        (
            "Mw_i",
            dimMass / dimMoles,
            gasMw[i]
        );

        invMwMix += Y[i] / Mw_i;
    }

    // Upper bound of the mixture molecular weight: 500 kg/kmol
    invMwMix = max
    (
        invMwMix,
        dimensionedScalar
        (
            "invMwMin",
            dimMoles / dimMass,
            scalar(1.0) / scalar(500.0)
        )
    );

    const volScalarField MwMix
    (
        IOobject
        (
            "MwMix", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        scalar(1.0) / invMwMix
    );

    // Density [kg/m3] from the ideal gas law at the reference pressure
    // (low Mach number), with the temperature bounded below by 298 K
    const dimensionedScalar R_univ
    (
        "R_univ",
        dimEnergy / dimMoles / dimTemperature,
        scalar(8314.46)   // [J/(kmol.K)]
    );

    const volScalarField T_bounded
    (
        IOobject
        (
            "T_bounded", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        max
        (
            T,
            dimensionedScalar("Tmin", dimTemperature, scalar(298.0))
        )
    );

    rho = (pRef_abs * MwMix) / (R_univ * T_bounded);

    rho = max
    (
        rho,
        dimensionedScalar("rhoMin", dimDensity, scalar(0.01))
    );
    rho.correctBoundaryConditions();

    // Specific heat capacity [J/(kg.K)]: Cp_mix = sum(Y_i*Cp_i)
    Cp = dimensionedScalar("zero", Cp.dimensions(), scalar(0.0));

    forAll(Y, i)
    {
        const dimensionedScalar Cp_i
        (
            "Cp_i",
            Cp.dimensions(),
            gasCp[i]
        );

        Cp += Y[i] * Cp_i;
    }

    Cp = max
    (
        Cp,
        dimensionedScalar("CpMin", Cp.dimensions(), scalar(100.0))
    );
    Cp.correctBoundaryConditions();

    // Thermal conductivity [W/(m.K)]: kappa_mix = sum(Y_i*kappa_i)
    kappa = dimensionedScalar("zero", kappa.dimensions(), scalar(0.0));

    forAll(Y, i)
    {
        const dimensionedScalar kappa_i
        (
            "kappa_i",
            kappa.dimensions(),
            gasKappa[i]
        );

        kappa += Y[i] * kappa_i;
    }

    kappa = max
    (
        kappa,
        dimensionedScalar("kappaMin", kappa.dimensions(), scalar(1e-6))
    );
    kappa.correctBoundaryConditions();
}
