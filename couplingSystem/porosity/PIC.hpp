#ifndef __PIC_hpp__ 
#define __PIC_hpp__

#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "porosity.hpp"


namespace pFlow::coupling
{


class PIC
: 
	public porosity
{
protected:

	

public:

	// type info
	TypeInfo("PIC");

	PIC(
		Foam::dictionary 		dict, 
		couplingMesh& 			cMesh, 
		MPI::centerMassField& 	centerMass, 
		MPI::realProcCMField& 	parDiam);

	virtual ~PIC() = default;

	add_vCtor
	(
		porosity,
		PIC,
		dictionary
	);

	bool calculatePorosity(Foam::volScalarField& alpha) override;

}; 

} // pFlow::coupling


#endif
