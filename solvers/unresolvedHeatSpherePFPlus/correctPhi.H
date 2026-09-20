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
 * @file correctPhi.H
 * @brief Evaluates and corrects the volumetric fluxes at the cell faces.
 */

// Calculate absolute flux
// from the mapped surface velocity
alphaRhoPhi = mesh.Sf() & Uf();

correctUphiBCs(U, alphaRhoPhi);

CorrectPhi
(
    U,
    alphaRhoPhi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(alphaRhoPhi, U);



