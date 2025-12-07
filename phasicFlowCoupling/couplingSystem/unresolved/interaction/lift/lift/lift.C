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

#include "fvc.H"

#include "lift.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"



namespace pFlow::coupling
{

lift::lift(const unresolvedCouplingSystem &uCS, const porosity &prsty)
:
    porosity_(prsty),
    printLift_(this->dict().getOrDefault("printLift", Foam::Switch(true)))
{
    surfaceRotationTorque_ = surfaceRotationTorque::create(uCS, prsty);
}

void lift::calculateLiftForceTorque
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleForce,
    Plus::realx3ProcCMField&        particleTorque
)
{
    calculateLiftForce(U, parVel, parRotVel, diameter, particleForce, particleTorque);

    surfaceRotationTorque_->calculateSurfaceTorque(U, parVel, parRotVel, diameter, particleTorque);
}


Foam::tmp<Foam::volVectorField> pFlow::coupling::lift::fluidVorticity(const Foam::volVectorField &U)
{    
    return Foam::fvc::curl(U);
}

const Foam::dictionary& lift::dict()const
{
    return lift::getDict(porosity_.uCS());
}


const Foam::dictionary& lift::getDict(const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("momentumInteraction").subDict("lift");
}

uniquePtr<lift> lift::create
(
    const unresolvedCouplingSystem& uCS, 
    const porosity& 				prsty
)
{
    const auto& liftDict = lift::getDict(uCS);

    //auto shapeName = uCS.shapeTypeName();
    auto liftType = liftDict.getOrDefault<Foam::word>("model", "none");
    
        
    if( couplingSystemvCtorSelector_.search(liftType))
    {
        Foam::Info<<"    Crearing lift force "<<Green_Text(liftType)<<" ...\n\n";
        return couplingSystemvCtorSelector_[liftType] (uCS, prsty);
    }
    else
    {
        if(Plus::processor::isMaster())
        {
            printKeys
            ( 
                fatalErrorInFunction << "Ctor Selector "<< liftType << " dose not exist"
                " for lift method in "<< liftDict.name()
                <<"\nAvaiable ones are: \n"
                ,
                couplingSystemvCtorSelector_
            )<<endl;
        }
        Plus::processor::abort(0);
    }

    return nullptr;
}

}