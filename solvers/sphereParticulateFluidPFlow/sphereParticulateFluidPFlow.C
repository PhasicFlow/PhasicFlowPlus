
// from OpenFoam
#include "fvCFD.H"
#include "fvModels.H"
#include "dynamicFvMesh.H"
#include "fvConstraints.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// from phasicFlow-coupling
#include "couplingSystem.hpp"



int main( int argc, char* argv[])
{
	#include "postProcess.H"

	#include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
	#include "createDyMControls.H"
    #include "createFields.H"

	
	
	return true;
}
