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

#include "foamPatchWall.hpp"
#include "streams.hpp"

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
	Foam::fileName caseName("./");


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

	Foam::Info<<names<<Foam::endl;

	if( auto pId = boundaries.findPatchID(patchName_); pId != -1 )
	{
		Foam::Info<<boundaries[patchName_]<<Foam::endl;
	}
	else
	{
		fatalExit;
	}

}