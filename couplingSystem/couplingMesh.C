
// from OpenFOAM



#include "couplingMesh.hpp"

void pFlow::coupling::couplingMesh::calculateBox()
{
	const auto& p = mesh_.points();
	auto lower = Foam::min(p);
	auto upper = Foam::max(p);

	meshBox_ = box(
		{lower[0], lower[1], lower[2]},
		{upper[0], upper[1], upper[2]});
}


pFlow::coupling::couplingMesh::couplingMesh(Foam::fvMesh& mesh) 
:
	mesh_(mesh),
	cellTreeSearch_(mesh.cellTree())
{
	calculateBox();
}



