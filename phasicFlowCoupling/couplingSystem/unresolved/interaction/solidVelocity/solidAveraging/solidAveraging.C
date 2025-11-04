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

template<typename DistributionType>
pFlow::coupling::solidAveraging<DistributionType>::solidAveraging
(
    const word& type,
    const couplingMesh& cMesh,
    Plus::realx3ProcCMField& particleField
)
:
    parAvField_(particleField)
{
}


template<typename DistributionType>
pFlow::uniquePtr<solidAveraging> pFlow::coupling::solidAveraging<DistributionType>::create
(
    const word& type,
    const word& distType,
    const couplingMesh& cMesh,
    Plus::realx3ProcCMField& particleField
)
{
    
    auto avType = angleBracketsNames2("solidAveraging", distType, type);

    if( wordvCtorSelector_.search(avType))
	{
		Foam::Info<<"    Crearing solid averaging "
                  << Green_Text(avType)
                  <<" for "<< Green_Text(particleField.name())<<"\n\n";
		return wordvCtorSelector_[avType] (type, distType, cMesh, particleField);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< avType << " dose not exist"
				" for drag method in "<< uDict.name()
				<<"\nAvaiable ones are: \n"
				,
				wordvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}

/*
    Us_
    (
        Foam::IOobject
	    (
	        "Us",
	        Foam::timeName(cMesh.mesh().time()),
	        cMesh.mesh(),
	        Foam::IOobject::NO_READ,
	        (distribute_?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
	    ),
    	cMesh.mesh(),
        Foam::dimensionedVector(Foam::dimVelocity, Foam::vector(0,0,0))
    )
*/