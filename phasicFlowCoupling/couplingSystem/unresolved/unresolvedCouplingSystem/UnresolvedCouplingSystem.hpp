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

#ifndef __UnresolvedCouplingSystem_hpp__
#define __UnresolvedCouplingSystem_hpp__


#include "unresolvedCouplingSystem.hpp"


namespace pFlow::coupling
{

template<typename DistributorType>
class UnresolvedCouplingSystem
:
	public unresolvedCouplingSystem
{
private:

	DistributorType 			cellDistribution_;

protected:

	DistributorType& cellDistribution()
	{
		return cellDistribution_;
	}

public:


	UnresolvedCouplingSystem(
		word shapeTypeName, 
		Foam::fvMesh& mesh,
		int argc, 
		char* argv[]);

	

	~UnresolvedCouplingSystem() override = default;


	const DistributorType& cellDistribution()const
	{
		return cellDistribution_;
	}

}; 

} // pFlow::coupling

#include "UnresolvedCouplingSystem.C"

#endif //__UnresolvedCouplingSystem_hpp__
