
// from OpenFOAM


#include "couplingMesh.hpp"
#include "processorPlus.hpp"
#include "procCommunicationPlus.hpp"
#include "streams.hpp"
#include "schedule.hpp"

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

    REPORT(1)<< Blue_Text("Bounding box has been updated.")<<END_REPORT;

}

void pFlow::coupling::couplingMesh::resetTree()const
{
    Foam::treeBoundBox bb(Foam::point(
                    meshBox_.minPoint().x(),
                    meshBox_.minPoint().y(),
                    meshBox_.minPoint().z()),
                Foam::point(
                    meshBox_.maxPoint().x(),
                    meshBox_.maxPoint().y(),
                    meshBox_.maxPoint().z()));

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
            treeBoundBoxExtend
            (
                bb,
                1.0e-3
            ),
            8,              // maxLevel
            10,             // leafsize
            6.0             // duplicity
        )
    );

    REPORT(1)<< Blue_Text("Search tree has been reset.")<<END_REPORT;
     
}

void pFlow::coupling::couplingMesh::mapParticles()
{
    const auto& cm = parCellIndex_.centerMass();
    const size_t numPar = cm.size();
    numInMesh_ = 0;
    
    #pragma ParallelRegion reduction(+:numInMesh_)
    for(size_t i = 0; i<numPar; i++)
    {
        auto cellId = findCellTree(cm[i], parCellIndex_[i]);
        parCellIndex_[i] = cellId;
        if( cellId >= 0 ) numInMesh_++;	
    }
}

pFlow::coupling::couplingMesh::couplingMesh
(
	const Foam::dictionary& dict,
    Foam::fvMesh& mesh, 
    const Plus::centerMassField& centerMass
) 
:
	mesh_(mesh),
    parCellIndex_
	(
		"parCellIndex",
		static_cast<Foam::label>(-1),
		centerMass,
		true
	)
{

    /// Mesh decomposition mode for locating points in the mesh
    /// Options are: 
    /// 	1) facePlanes
    /// 	2) cellTets
    /// 	3) faceDiagonalTriangles
    /// 	4) faceCenterTriangles
	Foam::word decompositionMode(
        lookupDict<Foam::word>(dict, "decompositionMode"));

    if(decompositionMode == "facePlanes") 
        cellDecompositionMode_ = Foam::polyMesh::FACE_PLANES;
    else if(decompositionMode == "cellTets")
        cellDecompositionMode_ = Foam::polyMesh::CELL_TETS;
    else if(decompositionMode == "faceDiagonalTriangles")
        cellDecompositionMode_ = Foam::polyMesh::FACE_DIAG_TRIS;
    else if(decompositionMode == "faceCenterTriangles")
        cellDecompositionMode_ = Foam::polyMesh::FACE_CENTRE_TRIS;
    else
    {
        fatalErrorInFunction<<
        "Wrong deompositionMode"<< decompositionMode <<
        " in dictionary "<< dict.name() <<endl;
        Plus::processor::abort(0);
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

    nCells_ = mesh_.nCells();
    calculateBox();
    resetTree();
}


void pFlow::coupling::couplingMesh::update()
{
    // for dynamic mesh, bounding box and search tree should be 
    // updated every time step
    if( dynamic() )
    {
        calculateBox();
        resetTree();
    }
    nCells_ = mesh_.nCells();
    mapParticles();
    reportNumInMesh();
}

void pFlow::coupling::couplingMesh::reportNumInMesh()const
{
    Plus::procCommunication proc;
	if( auto [numInMeshAll, success] = proc.collectAllToMaster(numInMesh_); success)
	{
		if(Plus::processor::isMaster())
		{
			int32 s=0;
			for(auto v:numInMeshAll) s += v;

			output<<Blue_Text("Particles located in processor meshes:") << 
			Yellow_Text(numInMeshAll)<<
			" => "<< Yellow_Text(s)<< endl;
		}
	}
}


bool pFlow::coupling::couplingMesh::pointInCell
(
	const Foam::point& p, 
	Foam::label celli
)const
{
	const Foam::labelList& f = mesh_.cells()[celli];
    const Foam::labelList& owner = mesh_.faceOwner();
    const Foam::vectorField& cf = mesh_.faceCentres();
    const Foam::vectorField& Sf = mesh_.faceAreas();
    

    Foam::scalar dist;
    forAll(f, facei)
    {
        Foam::label nFace = f[facei];
        Foam::vector proj = p - cf[nFace];
        
        if (owner[nFace] != celli)
        {
        	dist = Sf[nFace] & proj;
        }
        else
        {
        	dist = -(Sf[nFace] & proj);
        }
        
        if(dist < 0 )
        {
        	return false;
        }
    }

    return true;
}

bool pFlow::coupling::couplingMesh::pointSphereInCell
(
	const Foam::point& p, 
	Foam::scalar rad, 
	Foam::label celli, 
	bool& sphereInCell
) const
{

	const Foam::labelList& f = mesh_.cells()[celli];
    const Foam::labelList& owner = mesh_.faceOwner();
    const Foam::vectorField& cf = mesh_.faceCentres();
    const Foam::vectorField& Sf = mesh_.faceAreas();
    
    sphereInCell = true;
    

    Foam::scalar dist;
    forAll(f, facei)
    {
        Foam::label nFace = f[facei];
        Foam::vector proj = p - cf[nFace];
        Foam::vector normal = normalised(Sf[nFace]);
        
        if (owner[nFace] != celli)
        {
        	dist = normal & proj;
        }
        else
        {
        	dist = -(normal & proj);
        }
        
        if(dist < rad ) sphereInCell = false;
        if(dist < 0 )
        {
        	return false;
        }
    }

    return true;
}



bool pFlow::coupling::couplingMesh::pointSphereInCell
(
	const Foam::point& p, 
	Foam::scalar smallRad, 
	Foam::scalar largeRad,
	Foam::label celli, 
	bool& smallInCell, 
	bool& largeInCell
)const
{
	
	const Foam::labelList& f = mesh_.cells()[celli];
    const Foam::labelList& owner = mesh_.faceOwner();
    const Foam::vectorField& cf = mesh_.faceCentres();
    const Foam::vectorField& Sf = mesh_.faceAreas();
    
    smallInCell = true;
    largeInCell = true;

    Foam::scalar dist;
    forAll(f, facei)
    {
        Foam::label nFace = f[facei];
        Foam::vector proj = p - cf[nFace];
        Foam::vector normal = normalised(Sf[nFace]);
        
        if (owner[nFace] != celli)
        {
        	dist = normal & proj;
        }
        else
        {
        	dist = -(normal & proj);
        }
        
        if(dist < smallRad ) smallInCell = false;
        if(dist < largeRad ) largeInCell = false;
        if(dist < 0 )
        {
        	return false;
        }
    }

    return true;
}

Foam::label 
pFlow::coupling::couplingMesh::findPointInCellTree
(
	const Foam::point& p, 
	Foam::label cellId
)const
{
	if (cellId == -1 || cellId >= nCells_ )
    {
        
        return cellTreeSearch_().findInside(p);
    }
    else
    {
        if (pointInCell(p, cellId))
        {
            return cellId;
        }
        else
        {
            return cellTreeSearch_().findInside(p);
        }
    }
}

Foam::label 
pFlow::coupling::couplingMesh::findPointSphereInCellTree
(
	const Foam::point& p,
	Foam::scalar rad,
	Foam::label cellId,
	bool& sphereInCell
)const
{
	// first find the cellId if not known
	if (cellId == -1 || cellId >= nCells_ )
    {
        cellId = cellTreeSearch_().findInside(p);

        // point is not in mesh
        if(cellId ==-1)
        {
        	sphereInCell = false;
        	return cellId;
        }
    }
   
    if(pointSphereInCell(p, rad, cellId, sphereInCell)) return cellId;
    
    // the point may have been moved to another cell, find new cell
    cellId = cellTreeSearch_().findInside(p);
    if(cellId ==-1 || cellId >= nCells_)
    {
    	sphereInCell = false;
    	return cellId;
    }
    if(pointSphereInCell(p, rad, cellId, sphereInCell)) return cellId;
    return -1;
    
}

Foam::label 
pFlow::coupling::couplingMesh::findPointSphereInCellTree
(
	const Foam::point& p,
	Foam::scalar radSmall,
	Foam::scalar radLarge,
	Foam::label cellId,
	bool& smallInCell,
	bool& largeInCell
)const
{
	// first find the cellId if not known
	if (cellId == -1 || cellId >= nCells_)
    {
        cellId = cellTreeSearch_().findInside(p);

        // point is not in mesh
        if(cellId ==-1)
        {
        	smallInCell = false;
        	largeInCell = false;
        	return cellId;
        }
    }
   
    if(pointSphereInCell(p, radSmall, radLarge, cellId, smallInCell, largeInCell)) return cellId;
    
    // the point may have been moved to another cell, find new cell
    cellId = cellTreeSearch_().findInside(p);
    if(cellId ==-1 || cellId >= nCells_)
    {
    	smallInCell = false;
        largeInCell = false;
        return cellId;
    }
    if(pointSphereInCell(p, radSmall, radLarge, cellId, smallInCell, largeInCell)) return cellId;
    return -1;
}


Foam::label
pFlow::coupling::couplingMesh::findCellTree
(
    const realx3& p, 
    Foam::label cellId
)const
{
    Foam::point pp (p.x(), p.y(), p.z());
    if (cellId == -1 || cellId >= nCells_)
    {
        return cellTreeSearch_().findInside(pp);
    }
    else
    {
        if (pointInCell(
            pp,
            cellId))
        {
            return cellId;
        }
        else
        {
            return cellTreeSearch_().findInside(pp);     
        }
    }
}

Foam::labelList pFlow::coupling::couplingMesh::findSphere
(
    Foam::label cellId, 
    Foam::scalar radius
)const
{
    const auto& targetCellCentre = mesh_.cellCentres()[cellId];

     Foam::treeBoundBox searchBox
    (
        targetCellCentre - Foam::vector(radius, radius, radius),
        targetCellCentre + Foam::vector(radius, radius, radius)
    );
    
    return  cellTreeSearch_().findBox(searchBox);
}

