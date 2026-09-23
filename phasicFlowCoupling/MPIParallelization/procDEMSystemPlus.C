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

#include "procDEMSystemPlus.hpp"
#include "procCommunicationPlus.hpp"
#include "dictionary.hpp"

pFlow::Plus::procDEMSystem::procDEMSystem
(
	word demSystemName,
	int argc, 
	char* argv[],
	bool requireRVel
)
{
	procCommunication proc;

	if(Plus::processor::isMaster())
	{	
		realx3 domainMin(0, 0, 0);
		realx3 domainMax(1, 1, 1);

		// Standard hierarchy of search directories for domainDict
		std::vector<word> possibleDirs = {
			"settings",
			"caseSetup",
			"constant",
			"system",
			"."
		};

		word foundDir = "";

		for (const auto& dir : possibleDirs)
		{
			fileSystem candidate(dir, "domainDict");

			if (candidate.exist())
			{
				dictionary domDict("domainDict", candidate);
				const dictionary& globalBox = domDict.subDict("globalBox");

				domainMin = globalBox.getVal<realx3>("min");
				domainMax = globalBox.getVal<realx3>("max");

				foundDir = dir;
				break;
			}
		}

		if (!foundDir.empty())
		{
			output << "\n[PhasicFlow Plus] Read domain boundaries from '" 
				   << foundDir << "/domainDict':\n"
				   << "    Min: (" << domainMin.x() << " " 
				   << domainMin.y() << " " << domainMin.z() << ")\n"
				   << "    Max: (" << domainMax.x() << " " 
				   << domainMax.y() << " " << domainMax.z() << ")\n" 
				   << endl;
		}
		else
		{
			fatalErrorInFunction
				<< "CRITICAL ERROR: Cannot find 'domainDict' file in any "
				<< "standard folder.\n"
				<< "Checked: settings/, caseSetup/, constant/, system/, "
				<< "and root.\n"
				<< "This file is mandatory for correct NBS tracking of "
				<< "particles.\n"
				<< "Simulation aborted to prevent silent geometry corruption."
				<< endl;
			proc.abort(0);
		}

		demSystem_ = DEMSystem::create(
			demSystemName, 
			procVector<box>(
				box(domainMin, domainMax)), 
			argc, 
			argv,
			requireRVel);	
	}

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
