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

#include "polydisperseDragBase.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::polydisperseDragBase::polydisperseDragBase
(
    const unresolvedCouplingSystem& uCS, 
    const porosity& 				prsty
)
:
    drag(uCS, prsty),
    diameterClass_("diameterClass", uCS.centerMass()),
    writeAverageDiameter_
    (
        this->dict().template lookupOrDefault<Foam::Switch>
        (
            "writeAverageDiameter",
            Foam::Switch(false)
        )
    ),
    averageDiameter_
    (
        Foam::IOobject
        (
            "averageDiameter",
            mesh().time().timeName(),
            mesh(),
            Foam::IOobject::NO_READ,
            (writeAverageDiameter_?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
        ),
        mesh(),
        Foam::scalar(0)
    )
{}
