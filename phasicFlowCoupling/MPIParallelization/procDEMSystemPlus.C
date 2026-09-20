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
#include "procVectorPlus.hpp"
#include "procCommunicationPlus.hpp"

// Standard C++ inclusions needed for reading the domainDict file below
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

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

		// Standard hierarchy of search paths for domain bounds
		std::vector<std::string> possiblePaths = {
			"settings/domainDict",
			"caseSetup/domainDict",
			"constant/domainDict",
			"system/domainDict",
			"domainDict"
		};

		std::ifstream dictFile;
		std::string foundPath = "";

		for (const auto& path : possiblePaths)
		{
			dictFile.open(path);
			if (dictFile.is_open())
			{
				foundPath = path;
				break;
			}
		}

		if (dictFile.is_open())
		{
			std::string line;
			while (std::getline(dictFile, line))
			{
				if (line.find("min") != std::string::npos && 
					line.find("(")   != std::string::npos)
				{
					size_t start = line.find("(");
					size_t end   = line.find(")");

					if (start != std::string::npos && end != std::string::npos)
					{
						std::istringstream iss(
							line.substr(start + 1, end - start - 1));
						real x, y, z;
						if (iss >> x >> y >> z) 
						{
							domainMin = realx3(x, y, z);
						}
					}
				}
				else if (line.find("max") != std::string::npos && 
						 line.find("(")   != std::string::npos)
				{
					size_t start = line.find("(");
					size_t end   = line.find(")");

					if (start != std::string::npos && end != std::string::npos)
					{
						std::istringstream iss(
							line.substr(start + 1, end - start - 1));
						real x, y, z;
						if (iss >> x >> y >> z) 
						{
							domainMax = realx3(x, y, z);
						}
					}
				}
			}
			dictFile.close();

			std::cout << "\n[PhasicFlow Plus] Read domain boundaries from '" 
					  << foundPath << "':\n"
					  << "    Min: (" << domainMin.x() << " " 
					  << domainMin.y() << " " << domainMin.z() << ")\n"
					  << "    Max: (" << domainMax.x() << " " 
					  << domainMax.y() << " " << domainMax.z() << ")\n" 
					  << std::endl;
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
