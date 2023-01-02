
// from OpenFOAM
#include "fvCFD.H"


#include "PIC.hpp"



pFlow::coupling::PIC::PIC(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam)
{

}


bool pFlow::coupling::PIC::calculatePorosity(Foam::volScalarField& alpha)
{
	
	auto solidVol = Foam::scalarField(
		alpha.size(), 
		static_cast<Foam::scalar>(0));


	for(size_t i=0; i<centerMass_.size(); i++)
	{
		
		auto cellId = cMesh_.findCell(centerMass_[i], parCellIndex_[i]);
		if( cellId >= 0 )
		{
			solidVol[cellId] += 
				static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));
		}
		parCellIndex_[i] = cellId;
	}


	alpha.field() = Foam::max(
		1 - solidVol/alpha.mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()) );

	//alpha.ref() = 1 - solidVol/alpha.mesh().V();

	// we may need to update boundary condditions
	// we also need to check if the old time step value is stored or not.

	return true;
}