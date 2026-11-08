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
