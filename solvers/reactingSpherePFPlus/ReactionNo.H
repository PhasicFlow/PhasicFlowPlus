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

// Time step limit from the reaction rate. The fractional change of the gas
// mass of a cell per time step,
//
//     Rx = max(|Sm|/(alpha*rho)) * deltaT,
//
// may not exceed maxReactionNo (CLC dictionary of system/fvSolution). The
// limit is applied after setDeltaT.H and only when adjustTimeStep is on.

{
    const scalarField rateField
    (
        mag(Sm.primitiveField())
      / max
        (
            coupling.alpha().primitiveField() * rho.primitiveField(),
            scalar(SMALL)
        )
    );

    const scalar reactionNoNum = gMax(rateField) * runTime.deltaTValue();

    Info << "Reaction No (max|Sm|/(alpha*rho) * dt) : "
         << reactionNoNum << endl;

    if
    (
        adjustTimeStep
     && reactionNoNum > maxReactionNo
     && reactionNoNum > SMALL
    )
    {
        const scalar reactionFactor = maxReactionNo / reactionNoNum;

        const scalar newDeltaT =
            min(runTime.deltaTValue() * reactionFactor, maxDeltaT);

        Info << "  [CLC] Reaction-driven deltaT reduction: "
             << runTime.deltaTValue() << " -> " << newDeltaT
             << "  (Rx = " << reactionNoNum << ",  limit = "
             << maxReactionNo << ")" << endl;

        runTime.setDeltaT(newDeltaT);
    }
}

// Time step limit from the consumption of each gas species. The fractional
// change of the mass fraction of a consumed species of a cell per time step,
//
//     Rx_i = max(Sp_i/(alpha*Mw_i)) * deltaT,
//
// may not exceed maxSpeciesReactionNo (CLC dictionary of system/fvSolution).
// Sp_i/Mw_i is the sink coefficient Sp_c/rho of YEqn.H. The limit is applied
// after the one above and only when adjustTimeStep is on.

{
    const PtrList<volScalarField>& speciesSuFields = coupling.speciesSu();
    const PtrList<volScalarField>& speciesSpFields = coupling.speciesSp();

    scalar speciesReactionNoNum = 0;

    forAll(gasSpeciesNames, i)
    {
        forAll(speciesSuFields, k)
        {
            if (speciesSuFields[k].name() == "Su_" + gasSpeciesNames[i])
            {
                const scalarField rateField
                (
                    speciesSpFields[k].primitiveField()
                  / max
                    (
                        coupling.alpha().primitiveField() * gasMw[i],
                        scalar(SMALL)
                    )
                );

                speciesReactionNoNum = max
                (
                    speciesReactionNoNum,
                    gMax(rateField) * runTime.deltaTValue()
                );
            }
        }
    }

    Info << "Species Reaction No (max Sp/(alpha*Mw) * dt) : "
         << speciesReactionNoNum << endl;

    if
    (
        adjustTimeStep
     && speciesReactionNoNum > maxSpeciesReactionNo
     && speciesReactionNoNum > SMALL
    )
    {
        const scalar reactionFactor =
            maxSpeciesReactionNo / speciesReactionNoNum;

        const scalar newDeltaT =
            min(runTime.deltaTValue() * reactionFactor, maxDeltaT);

        Info << "  [CLC] Species-driven deltaT reduction: "
             << runTime.deltaTValue() << " -> " << newDeltaT
             << "  (Rx_i = " << speciesReactionNoNum << ",  limit = "
             << maxSpeciesReactionNo << ")" << endl;

        runTime.setDeltaT(newDeltaT);
    }
}
