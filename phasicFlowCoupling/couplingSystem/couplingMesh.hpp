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

#ifndef __couplingMesh_hpp__ 
#define __couplingMesh_hpp__

// from OpneFOAM
#include "fvMesh.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// phasicFlow
#include "box.hpp"


namespace pFlow::coupling
{

class couplingMesh
{
protected:

	const Foam::fvMesh& 					mesh_;

	Foam::polyMesh::cellDecomposition 		cellDecompMode_;

	mutable uniquePtr<Foam::indexedOctree<Foam::treeDataCell>>
		cellTreeSearch_ = nullptr;

	box 			meshBox_;
	
	void calculateBox();

	Foam::label findCellSeed(
		const Foam::point& loc,
    	const Foam::label seedCellId);

public:

	couplingMesh(
		Foam::fvMesh& mesh, 
		const Foam::polyMesh::cellDecomposition 
				decompMode = Foam::polyMesh::FACE_PLANES);

	~couplingMesh()=default;

	inline
	auto& mesh()
	{
		return mesh_;
	}

	inline
	const auto& mesh()const
	{
		return mesh_;
	}

	inline
	const auto& meshBox()const
	{
		return meshBox_;
	}

	
	Foam::label 
	findCell(const realx3& p, Foam::label cellCheck);

	Foam::label
	findCellTree(const realx3& p, Foam::label cellId);




};

}



#endif //__couplingMesh_hpp__
