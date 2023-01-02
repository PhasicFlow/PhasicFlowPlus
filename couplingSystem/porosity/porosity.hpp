#ifndef __porosity_hpp__ 
#define __porosity_hpp__

// from OpenFOAM
#include "dictionary.H"

// from phasicFlow
#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "couplingMesh.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{


class porosity
{
protected:

	real 							alphaMin_;

	couplingMesh& 					cMesh_;

	MPI::centerMassField& 			centerMass_;

	MPI::realProcCMField&  			particleDiameter_;

	MPI::procCMField<Foam::label> 	parCellIndex_;

public:

	// type info
	TypeInfo("porosity");

	porosity(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	virtual ~porosity() = default;

	create_vCtor
	(
		porosity,
		dictionary,
		(
			Foam::dictionary		dict, 
			couplingMesh& 			cMesh, 
			MPI::centerMassField& 	centerMass, 
			MPI::realProcCMField& 	parDiam
		),
		(dict, cMesh, centerMass, parDiam)
	);

	inline
	real alphaMin()const
	{
		return alphaMin_;
	}

	virtual
	bool calculatePorosity(Foam::volScalarField& alpha) = 0;

	static
	uniquePtr<porosity> create(
		Foam::dictionary		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

}; 

} // pFlow::coupling


#endif
