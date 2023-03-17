
// from OpenFOAM



#include "couplingMesh.hpp"

void pFlow::coupling::couplingMesh::calculateBox()
{
	const auto& p = mesh_.points();
	auto lower = Foam::min(p);
	auto upper = Foam::max(p);

	meshBox_ = box(
		{static_cast<real>(lower[0]), 
         static_cast<real>(lower[1]), 
         static_cast<real>(lower[2])},
		{static_cast<real>(upper[0]), 
         static_cast<real>(upper[1]), 
         static_cast<real>(upper[2])});
}

Foam::label 
pFlow::coupling::couplingMesh::findCell
(
	const realx3& p, 
	Foam::label cellCheck
)
{
	return mesh_.findCell(Foam::point(p.x(), p.y(), p.z()));
}


Foam::label 
pFlow::coupling::couplingMesh::findCellSeed
(
	const Foam::point& loc,
    const Foam::label seedCellId
)
{
	if (mesh_.pointInCell(loc, seedCellId, cellDecompMode_))
    {
        return seedCellId;
    }

    Foam::label  curCelli = seedCellId;
    Foam::scalar nearestDistSqr = magSqr(mesh_.cellCentres()[seedCellId] - loc);

    while(true)
    {
        // Try neighbours of curCelli

        const auto& cFaces = mesh_.cells()[curCelli];

        Foam::label nearestCelli = -1;

        forAll(cFaces, i)
        {
            Foam::label facei = cFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                Foam::label celli = mesh_.faceOwner()[facei];
                if (celli == curCelli)
                {
                    celli = mesh_.faceNeighbour()[facei];
                }

                // Check if this is the correct cell
                if (mesh_.pointInCell(loc, celli, cellDecompMode_))
                {
                    return celli;
                }

                // Also calculate the nearest cell
                Foam::scalar distSqr = Foam::magSqr(mesh_.cellCentres()[celli] - loc);

                if (distSqr < nearestDistSqr)
                {
                    nearestDistSqr = distSqr;
                    nearestCelli = celli;
                }
            }
        }

        if (nearestCelli == -1)
        {
            return -1;
        }

        // Continue with the nearest cell
        curCelli = nearestCelli;
    }

    return -1;
}


pFlow::coupling::couplingMesh::couplingMesh(
	Foam::fvMesh& mesh,
	const Foam::polyMesh::cellDecomposition decompMode) 
:
	mesh_(mesh),
	cellDecompMode_(decompMode)
{

	if
    (
        decompMode == Foam::polyMesh::FACE_DIAG_TRIS
     || decompMode == Foam::polyMesh::CELL_TETS)
    {
        (void)mesh.tetBasePtIs();
    }

    if(!cellTreeSearch_)
    {
    	cellTreeSearch_.reset
        (
            new Foam::indexedOctree<Foam::treeDataCell>
            (
                Foam::treeDataCell
                (
                    false,      // not cache bb
                    mesh,
                    decompMode   // use tet-decomposition for any inside test
                ),
                Foam::treeBoundBox(mesh.points()).extend(1e-4),
                8,              // maxLevel
                10,             // leafsize
                6.0             // duplicity
            )
        );
    }

	calculateBox();
}

Foam::label
pFlow::coupling::couplingMesh::findCellTree(
	const realx3& p, 
	Foam::label cellId)
{
	if (cellId == -1)
    {
        return cellTreeSearch_().findInside(
        	Foam::point(p.x(), p.y(), p.z())
        	);
    }
    else
    {
        return findCellSeed(
        	Foam::point(p.x(), p.y(), p.z()), cellId
        	);
    }
}

