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

// Gas species transport in kinematic form (all terms divided by rho).
//
// Su is the net explicit source of a species [kg/(m3.s)] and Sp the implicit
// sink coefficient [kg/(kmol.s)] of the DEM reactions, both from the previous
// DEM step and frozen over the PIMPLE iterations. With Sp_c = Sp*rho/Mw and
// Yold = Y.oldTime(), each equation reads
//
//     ddt(alpha*Y) + div(alphaRhoPhi*Y) - laplacian(alpha*Deff, Y)
//         = (Su + Sp_c*Yold)/rho - Sp_c/rho*Y - (Sm/rho - Sv)*Y
//
// (Picard linearisation of the sink). The last term is the difference
// between the mass source Sm/rho and the volume source Sv of the continuity
// equation. All species except the inert one are transported; the inert
// species closes the sum to one.

{
    // Effective mass diffusivity [m2/s]
    const volScalarField Deff
    (
        "Deff",
        turbulence->nu() / Sc + turbulence->nut() / Sct
    );

    const volScalarField& alpha = coupling.alpha();

    // The lists follow couplingProperties/productGas, which need not match
    // the order of Y, so each species is looked up by the name
    // "Su_<species>" below.
    const PtrList<volScalarField>& speciesSuFields = coupling.speciesSu();
    const PtrList<volScalarField>& speciesSpFields = coupling.speciesSp();

    // Sources of a species that takes part in no reaction
    const volScalarField zeroSu
    (
        IOobject
        (
            "zeroSu", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroSu", dimDensity/dimTime, scalar(0))
    );

    const volScalarField zeroSp
    (
        IOobject
        (
            "zeroSp", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroSp", dimMass/dimMoles/dimTime, scalar(0))
    );

    // Sum of the transported fractions
    volScalarField Yt
    (
        IOobject
        (
            "Yt", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Yt", dimless, scalar(0))
    );

    forAll(Y, i)
    {
        if (i == inertIndex) continue;

        volScalarField& Yi = Y[i];

        // Zero sources when the species is not in productGas
        label rxnIdx = -1;
        forAll(speciesSuFields, k)
        {
            if (speciesSuFields[k].name() == "Su_" + gasSpeciesNames[i])
            {
                rxnIdx = k;
                break;
            }
        }

        const volScalarField& SuYi =
            (rxnIdx >= 0) ? speciesSuFields[rxnIdx] : zeroSu;
        const volScalarField& SpYi =
            (rxnIdx >= 0) ? speciesSpFields[rxnIdx] : zeroSp;

        // Sp_c = Sp*rho/Mw [kg/(m3.s)]
        const dimensionedScalar Mw_i
        (
            "Mw_" + gasSpeciesNames[i],
            dimMass / dimMoles,
            gasMw[i]
        );

        const volScalarField Sp_Corrected
        (
            "Sp_" + gasSpeciesNames[i] + "_corr",
            SpYi * rho / Mw_i
        );

        const volScalarField Su_Effective
        (
            "Su_" + gasSpeciesNames[i] + "_eff",
            SuYi + Sp_Corrected * Yi.oldTime()
        );

        fvScalarMatrix YiEqn
        (
            fvm::ddt(alpha, Yi)
          + fvm::div(alphaRhoPhi, Yi)
          - fvm::laplacian(alpha * Deff, Yi)
          ==
            Su_Effective / rho
          - fvm::Sp(Sp_Corrected / rho, Yi)

            // Mass source Sm/rho against the volume source Sv of continuity
          - fvm::SuSp(Sm / rho - Sv, Yi)

          + fvOptions(Yi)
        );

        YiEqn.relax();
        fvOptions.constrain(YiEqn);
        YiEqn.solve(mesh.solver(Yi.name()));
        fvOptions.correct(Yi);

        Yi.max(scalar(0));
        Yi.min(scalar(1));
        Yi.correctBoundaryConditions();

        Yt += Yi;
    }

    // Scale factor min(1, 1/Yt) keeps the sum of the transported fractions
    // at or below one, so that the inert species stays non-negative
    volScalarField Yscale
    (
        "Yscale",
        min
        (
            scalar(1.0),
            scalar(1.0) / (Yt + dimensionedScalar("eps", dimless, 1e-15))
        )
    );

    forAll(Y, i)
    {
        if (i == inertIndex) continue;

        Y[i] *= Yscale;
        Y[i].max(scalar(0));
        Y[i].min(scalar(1));
        Y[i].correctBoundaryConditions();
    }

    Y[inertIndex] = scalar(1) - (Yt * Yscale);
    Y[inertIndex].max(scalar(0));
    Y[inertIndex].min(scalar(1));
    Y[inertIndex].correctBoundaryConditions();

    // Yscale < 1 means that the transported fractions summed to more than
    // one before the scaling
    if (pimple.finalIter())
    {
        const scalar minYscale = min(Yscale).value();

        if (minYscale < scalar(1.0) - SMALL)
        {
            ++YscaleActivations;

            Info << "  [YEqn INFO] Active-species normalisation engaged "
                 << "this step (min Yscale = " << minYscale
                 << "). Cumulative occurrences since run start: "
                 << YscaleActivations
                 << endl;
        }
    }

    // Sum of all species fractions
    if (pimple.finalIter())
    {
        volScalarField Ytotal("Ytotal", Y[0]);
        for (label g = 1; g < Y.size(); ++g) Ytotal += Y[g];

        const scalar Ymin = min(Ytotal).value();
        const scalar Ymax = max(Ytotal).value();

        if (Ymax > scalar(1.02) || Ymin < scalar(-0.02))
        {
            Info << "  [YEqn WARNING] Species sum out of bounds: "
                 << "min(sum Y) = " << Ymin
                 << ",  max(sum Y) = " << Ymax << nl
                 << "  -> Check preExpFactor values and stoichiometric "
                 << "coefficients."
                 << endl;
        }
    }
}
