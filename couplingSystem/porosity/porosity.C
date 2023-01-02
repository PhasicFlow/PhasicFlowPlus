

#include "porosity.hpp"
#include "processor.hpp"

pFlow::coupling::porosity::porosity(
	Foam::dictionary		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	alphaMin_(dict.lookup<real>("alphaMin")),
	cMesh_(cMesh),
	centerMass_(centerMass),
	particleDiameter_(parDiam),
	parCellIndex_(
		"parCellIndex",
		static_cast<Foam::label>(-1),
		centerMass,
		true)
{

}


pFlow::uniquePtr<pFlow::coupling::porosity> 
	pFlow::coupling::porosity::create(
		Foam::dictionary		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam)
{
	auto method = dict.lookup<Foam::word>("method");
	if( dictionaryvCtorSelector_.search(method))
	{
		return dictionaryvCtorSelector_[method] (dict, cMesh, centerMass, parDiam);
	}
	else
	{
		if(MPI::processor::isMaster())
		{
			printKeys
			( 
				fatalErrorInFunction << "Ctor Selector "<< method << " dose not exist"
				" for porosity method in "<< dict.name()
				<<"\nAvaiable ones are: \n"
				,
				dictionaryvCtorSelector_
			)<<endl;
		}
		MPI::processor::abort(0);
	}

	return nullptr;
}