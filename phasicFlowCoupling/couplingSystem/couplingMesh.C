
// from OpenFOAM


#include "couplingMesh.hpp"
#include "processor.hpp"
#include "streams.hpp"

void pFlow::coupling::couplingMesh::calculateBox()const
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

    Foam::Info<< blueText("Bounding box has been updated.")<<Foam::endl;

}

void pFlow::coupling::couplingMesh::resetTree()const
{
    cellTreeSearch_.reset
    (
        new Foam::indexedOctree<Foam::treeDataCell>
        (
            Foam::treeDataCell
            (
                true,      // not cache bb
                mesh_,
                cellDecompositionMode_   // use tet-decomposition for any inside test
            ),
            Foam::treeBoundBox
            (
                Foam::point(
                    meshBox_.minPoint().x(),
                    meshBox_.minPoint().y(),
                    meshBox_.minPoint().z()),
                Foam::point(
                    meshBox_.maxPoint().x(),
                    meshBox_.maxPoint().y(),
                    meshBox_.maxPoint().z())
            ).extend(1e-3),
            8,              // maxLevel
            10,             // leafsize
            6.0             // duplicity
        )
    );

    Foam::Info<< blueText("Search tree has been reset.")<<Foam::endl;
}


pFlow::coupling::couplingMesh::couplingMesh
(
	const Foam::dictionary& dict,
    Foam::fvMesh& mesh
) 
:
	mesh_(mesh),
    meshStatic_(!mesh.dynamic()),
    domainExpansionRatio_
    (
        Foam::max(dict.lookup<Foam::scalar>("domainExpansionRatio"), 0.5)
    ),
    domainUpdateInterval_
    (
        dict.lookup<Foam::scalar>("domainUpdateInterval")
    ),
    decompositionMode_
    (
        dict.lookup<Foam::word>("decompositionMode")
    )
{

    if(decompositionMode_ == "facePlanes") 
        cellDecompositionMode_ = Foam::polyMesh::FACE_PLANES;
    else if(decompositionMode_ == "cellTets")
        cellDecompositionMode_ = Foam::polyMesh::CELL_TETS;
    else if(decompositionMode_ == "faceDiagonalTriangles")
        cellDecompositionMode_ = Foam::polyMesh::FACE_DIAG_TRIS;
    else if(decompositionMode_ == "faceCenterTriangles")
        cellDecompositionMode_ = Foam::polyMesh::FACE_CENTRE_TRIS;
    else
    {
        fatalErrorInFunction<<
        "Wrong deompositionMode"<< decompositionMode_ <<
        " in dictionary "<< dict.name() <<endl;
        MPI::processor::abort(0);
        return;
    }

	if
    (
        cellDecompositionMode_ == Foam::polyMesh::FACE_DIAG_TRIS
     || cellDecompositionMode_ == Foam::polyMesh::CELL_TETS
    )
    {
        (void)mesh.tetBasePtIs();
    }
    resetTree();
    calculateBox();
}


void pFlow::coupling::couplingMesh::update(Foam::scalar t, Foam::scalar fluidDt)
{
    
    checkForDomainUpdate(t, fluidDt);
    
    // for the first time, they should be constructed anyway
    if(!firstConstruction_) 
    {
        firstConstruction_ = true;
        calculateBox();
        resetTree();
    }

    // for dynamic mesh, bounding box and search tree should be 
    // updated every time step
    if( dynamic() )
    {
        calculateBox();
        resetTree();
    }
    
}

bool pFlow::coupling::couplingMesh::checkForDomainUpdate
(
    Foam::scalar t, 
    Foam::scalar fluidDt,
    bool insideFluidLoop
)
{
    if(insideFluidLoop) t -= fluidDt;

    if( !firstConstruction_ )
    {
        lastTimeUpdated_ = t;
        return true;
    }

    if( abs(t-lastTimeUpdated_) < 0.98*fluidDt )
    {
        lastTimeUpdated_ = t;
        return true;
    }
    
    if( abs(t-(lastTimeUpdated_+domainUpdateInterval_)) < 0.98*fluidDt)
    {
        lastTimeUpdated_ = t;
        return true;
    }

    return false;
}


Foam::label
pFlow::coupling::couplingMesh::findCellTree
(
    const realx3& p, 
    Foam::label cellId
)const
{
    if (cellId == -1)
    {
        //output<<"cellId \n";
        return cellTreeSearch_().findInside(
            Foam::point(p.x(), p.y(), p.z())
            );
    }
    else
    {
        if (mesh_.pointInCell(
            Foam::point(p.x(), p.y(), p.z()),
            cellId,
            cellDecompositionMode_))
        {
            return cellId;
        }
        else
        {
            return cellTreeSearch_().findInside(
            Foam::point(p.x(), p.y(), p.z())
            );     
        }
    }
    /*Foam::point pt(p.x(), p.y(), p.z());

    if(cellId != -1 && 
        mesh_.pointInCell(pt, cellId, cellDecompositionMode_)
        ) return cellId;

    return mesh_.findCell(pt, Foam::polyMesh::CELL_TETS);*/

}

/*Foam::label 
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
}*/