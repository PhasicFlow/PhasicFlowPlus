/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------  
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

// from OpenFOAM
#include "Time.H"
#include "polyMesh.H"
#include "meshTriangulation.H"

#include "foamPatchWall.hpp"
#include "streams.hpp"


void addFaceTriangles(
	const Foam::face& f, 
	const Foam::Field<Foam::point>& points, 
	std::vector<pFlow::realx3x3>& triangles)
{
	auto n = f.size();

	auto nTris = n-2;
	pFlow::realx3 p0,p1,p2;
	auto ofP0 = points[f[0]];
	
	p0 = {ofP0.x(), ofP0.y(), ofP0.z()};



	for(Foam::label i=0; i<nTris; i++)
	{
		auto ofP1 = points[f[i+1]];
		auto ofP2 = points[f[i+2]];

		p1 = {ofP1.x(), ofP1.y(), ofP1.z()};
		p2 = {ofP2.x(), ofP2.y(), ofP2.z()};
		triangles.push_back({p0,p1,p2});
	}
}

pFlow::coupling::foamPatchWall::foamPatchWall()
{
}

pFlow::coupling::foamPatchWall::foamPatchWall(
  const dictionary& dict)
:
  Wall(dict),
  patchName_(dict.getVal<word>("patch"))
{


	output<<"foamPatchWall is created"<<endl;

	output<< "Creating time for OpenFOAM mesh ...\n"<<endl;
	
	Foam::fileName root(CWD().wordPath());
	Foam::fileName caseName("");


	Foam::Time runTime(root, caseName);

	output<< "Creating mesh ...\n"<<endl;
	Foam::polyMesh mesh
	(
	    Foam::IOobject
	    (
	        Foam::polyMesh::defaultRegion,
	        runTime.timeName(),
	        runTime,
	        Foam::IOobject::MUST_READ
	    )
	);

	auto& boundaries = mesh.boundaryMesh();

	auto names = boundaries.names();

	Foam::label pId = boundaries.findPatchID(patchName_);
	if( pId == -1 )
	{
		fatalExit;
	}

	auto& patch = boundaries[pId];	
	auto& lfaces = patch.localFaces();
	auto& lPoints = patch.localPoints();

	forAll(lfaces, fi)
	{
		addFaceTriangles(lfaces[fi], lPoints, triangles_);
	}

	REPORT(1)<<"Number of triagnles in patch "<< greenText(patch.name())<<
	" is " << yellowText(triangles_.size())<<endl;
	
}
