/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
License
    This file is partially part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This file is partially part of phasicFlow. 

    It is a free software for simulating 
    granular and multiphase flows. You can redistribute it and/or modify it under
    the terms of GNU General Public License v3 or any other later versions. 
 
    phasicFlow is distributed to help others in their research in the field of 
    granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#include "alphaTurbulentTransportModels.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTurbulenceModelTypes
(
    volScalarField,
    geometricOneField,
    incompressibleTurbulenceModel,
    AlphaIncompressibleTurbulenceModel,
    transportModel
);

makeBaseTurbulenceModel
(
    volScalarField,
    geometricOneField,
    incompressibleTurbulenceModel,
    AlphaIncompressibleTurbulenceModel,
    transportModel
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Stokes.H"
makeLaminarModel(Stokes);


#include "generalizedNewtonian.H"
makeLaminarModel(generalizedNewtonian);

#include "Maxwell.H"
makeLaminarModel(Maxwell);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "SpalartAllmaras.H"
makeRASModel(SpalartAllmaras);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "RNGkEpsilon.H"
makeRASModel(RNGkEpsilon);

#include "realizableKE.H"
makeRASModel(realizableKE);

#include "LaunderSharmaKE.H"
makeRASModel(LaunderSharmaKE);

#include "kEpsilonPhitF.H"
makeRASModel(kEpsilonPhitF);

#include "kOmega.H"
makeRASModel(kOmega);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSAS.H"
makeRASModel(kOmegaSSTSAS);

#include "kOmegaSSTLM.H"
makeRASModel(kOmegaSSTLM);

#include "LRR.H"
makeRASModel(LRR);

#include "SSG.H"
makeRASModel(SSG);

#include "EBRSM.H"
makeRASModel(EBRSM);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "WALE.H"
makeLESModel(WALE);

#include "kEqn.H"
makeLESModel(kEqn);

#include "dynamicKEqn.H"
makeLESModel(dynamicKEqn);

#include "dynamicLagrangian.H"
makeLESModel(dynamicLagrangian);

#include "sigma.H"
makeLESModel(sigma);

#include "SpalartAllmarasDES.H"
makeLESModel(SpalartAllmarasDES);

#include "SpalartAllmarasDDES.H"
makeLESModel(SpalartAllmarasDDES);

#include "SpalartAllmarasIDDES.H"
makeLESModel(SpalartAllmarasIDDES);

#include "DeardorffDiffStress.H"
makeLESModel(DeardorffDiffStress);

#include "kOmegaSSTDES.H"
makeLESModel(kOmegaSSTDES);

#include "kOmegaSSTDDES.H"
makeLESModel(kOmegaSSTDDES);

#include "kOmegaSSTIDDES.H"
makeLESModel(kOmegaSSTIDDES);
 

// ************************************************************************* //
