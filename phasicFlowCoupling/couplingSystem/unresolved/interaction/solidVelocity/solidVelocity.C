#include "solidVelocity.hpp"
#include "couplingMesh.hpp"


pFlow::coupling::solidVelocity::solidVelocity
(
    const word &type, 
    const Plus::realx3ProcCMField &parVel, 
    const couplingMesh &cMesh
)
:
    Up_(parVel),
    distribute_(type == "cell"),
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
{
}
