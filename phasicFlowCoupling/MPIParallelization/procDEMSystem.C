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

#include "procDEMSystem.hpp"
#include "procVector.hpp"
#include "procCommunication.hpp"

pFlow::MPI::procDEMSystem::procDEMSystem
(
	word demSystemName,
	int argc, 
	char* argv[]
)
{
	if(MPI::processor::isMaster())
	{	
		demSystem_ = DEMSystem::create(
			demSystemName, 
			procVector<box>(
				box
				(
					realx3(0), 
					realx3(1)
				)), 
			argc, 
			argv);	
	}

	procCommunication proc;
	
	real startT;
	if(demSystem_)
	{
		startT = demSystem_->Control().time().startTime();
	}
	else
	{
		startT = 0;
	}

	if(!proc.distributeMasterToAll(startT, startTime_))
	{
		fatalErrorInFunction<< "could not get start time"<<endl;
		proc.abort(0);
	}
}
