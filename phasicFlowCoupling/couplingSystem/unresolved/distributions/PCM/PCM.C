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

#include "PCM.hpp"

void pFlow::coupling::PCM::constructLists
(
    const Foam::scalar searchLen, 
    const Foam::label maxLayers
)
{
}

void pFlow::coupling::PCM::constructLists(const Foam::label maxLayers)
{
}

pFlow::coupling::PCM::PCM(
    const Foam::dictionary &parrentDict,
    const couplingMesh &cMesh,
    const Plus::centerMassField &centerMass)
: 
    distributionBase(false, parrentDict, cMesh, centerMass)
{}

pFlow::coupling::PCM::PCM
( 
    const couplingMesh &cMesh, 
    const Plus::centerMassField &centerMass
)
:
    distributionBase(false, cMesh, centerMass)
{
}

void pFlow::coupling::PCM::updateWeights(const Plus::procCMField<real> &parDiameter)
{
}

void pFlow::coupling::PCM::smoothenField(Foam::volVectorField &field) const
{
}

void pFlow::coupling::PCM::smoothenField(Foam::volScalarField &field) const
{
}
