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


template<typename DistributorType>
void pFlow::coupling::drag<DistributorType>::setSuSpToZero()
{
	// initialize all terms to zero
	Su_ = Foam::dimensionedVector(
    	Su_.name(),
    	Su_.dimensions(),
    	Foam::vector(0,0,0));

	Sp_ = Foam::dimensionedScalar(
		Sp_.name(), 
		Sp_.dimensions(), 
		0.0);
}


template<typename DistributorType>
pFlow::coupling::drag<DistributorType>::drag
(
	const UnresolvedCouplingSystem<DistributorType>& uCS, 
	const porosity& 				prsty
)
:
	porosity_(prsty),
	cellDistribution_(uCS.cellDistribution()),
	p_
	(
		porosity_.mesh().lookupObject<Foam::volScalarField>("p")
	),
	Su_
	(
		Foam::IOobject
		(
			"Su",
			porosity_.mesh().time().timeName(),
			porosity_.mesh(),
			Foam::IOobject::READ_IF_PRESENT,
			Foam::IOobject::AUTO_WRITE
		),
    	porosity_.mesh(),
    	Foam::dimensionedVector
    	(
	    	"Su",
	    	Foam::dimensionSet(1,-2,-2,0,0),
	    	Foam::vector(0,0,0)
	    )
	),
    Sp_
    (
    	Foam::IOobject
	    (
	        "Sp",
	        porosity_.mesh().time().timeName(),
	        porosity_.mesh(),
	        Foam::IOobject::READ_IF_PRESENT,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	porosity_.mesh(),
    	Foam::dimensionedScalar
    	(
    		"Sp", 
    		Foam::dimensionSet(1,-3,-1,0,0), 
    		0.0
    	)
    ),
    isCompressible_( p_.dimensions() == Foam::dimPressure)
{

}

template<typename DistributorType>
Foam::tmp<Foam::volVectorField> 
pFlow::coupling::drag<DistributorType>::pressureGradient(const Foam::volScalarField& rho)const
{
	if(isCompressible_)
		return Foam::fvc::grad(p_);
	else
		return Foam::fvc::grad(p_*rho);
}


template<typename DistributorType>
pFlow::uniquePtr<pFlow::coupling::drag<DistributorType>> 
pFlow::coupling::drag<DistributorType>::create
(
	const UnresolvedCouplingSystem<DistributorType>& uCS, 
	const porosity& 				prsty
)
{
	const auto& uDict = uCS.unresolvedDict();
	const auto& dragDict = uDict.subDict("drag");

	auto shapeName = uCS.shapeTypeName();
	auto useCellDist = dragDict.template lookup<Foam::Switch>("cellDistribution");
	auto distName = DistributorType::TYPENAME();
	auto dType = dragDict.template lookup<Foam::word>("type");
	
	Foam::word dragType;

	if(useCellDist)
	{
		dragType = angleBracketsNames3(shapeName+"Drag", distName, dType,"withDistribution");
	}
	else
	{
		dragType = angleBracketsNames3(shapeName+"Drag", distName, dType, "noDistribution");
	}

	
	if( couplingSystemvCtorSelector_.search(dragType))
	{
		Foam::Info<<"    Crearing drag force "<<Green_Text(dragType)<<" ...\n\n";
		return couplingSystemvCtorSelector_[dragType] (uCS, prsty);
	}
	else
	{
		if(Plus::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< dragType << " dose not exist"
				" for drag method in "<< uDict.name()
				<<"\nAvaiable ones are: \n"
				,
				couplingSystemvCtorSelector_
			)<<endl;
		}
		Plus::processor::abort(0);
	}

	return nullptr;
}