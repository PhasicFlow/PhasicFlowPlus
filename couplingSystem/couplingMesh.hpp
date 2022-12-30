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

#include "box.hpp"


namespace pFlow::coupling
{

class couplingMesh
{
protected:

	Foam::fvMesh& 							 mesh_;

	const 
	Foam::indexedOctree<Foam::treeDataCell>& cellTreeSearch_;

	box 			meshBox_;


	
	void calculateBox();

public:

	couplingMesh(Foam::fvMesh& mesh);

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

	/*inline
	const auto& foamTime()const
	{
		return mesh_.Time();
	}*/

	
	Foam::label 
	findCell(const realx3& p, Foam::label cellCheck)
	{
		#if useDouble
			const auto& sample = reinterpret_cast<const Foam::point&>(p);
		#else
			auto sample = Foam::point(p.x(), p.y(), p.z());
		#endif

		return -1;
	}



};

}



#endif //__couplingMesh_hpp__
