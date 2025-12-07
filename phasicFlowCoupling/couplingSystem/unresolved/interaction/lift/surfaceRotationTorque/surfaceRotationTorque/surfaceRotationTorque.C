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


#include "surfaceRotationTorque.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::surfaceRotationTorque::surfaceRotationTorque
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    porosity_(prsty)
{
}

const Foam::dictionary &pFlow::coupling::surfaceRotationTorque::dict() const
{
    return surfaceRotationTorque::getDict(porosity_.uCS());
}

const Foam::dictionary &pFlow::coupling::surfaceRotationTorque::getDict(const unresolvedCouplingSystem &uCS)
{
    return uCS.unresolvedDict().subDict("momentumInteraction").subDict("lift");
}


pFlow::uniquePtr<pFlow::coupling::surfaceRotationTorque> 
    pFlow::coupling::surfaceRotationTorque::create
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
{
    
    const auto& liftDict = getDict(uCS);

    //auto shapeName = uCS.shapeTypeName();
    auto surfaceRotationTorqueType = liftDict.getOrDefault<Foam::word>("surfaceRotationTorque", "none");
    
        
    if( couplingSystemvCtorSelector_.search(surfaceRotationTorqueType))
    {
        Foam::Info<<"    Crearing surfaceRotationTorque model "<<Green_Text(surfaceRotationTorqueType)<<" ...\n\n";
        return couplingSystemvCtorSelector_[surfaceRotationTorqueType] (uCS, prsty);
    }
    else
    {
        if(Plus::processor::isMaster())
        {
            printKeys
            ( 
                fatalErrorInFunction << "Ctor Selector "<< surfaceRotationTorqueType << " dose not exist"
                " for surfaceRotationTorque method in "<< liftDict.name()
                <<"\nAvaiable ones are: \n"
                ,
                couplingSystemvCtorSelector_
            )<<endl;
        }
        Plus::processor::abort(0);
    }
    return nullptr;
}