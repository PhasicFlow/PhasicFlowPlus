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

#ifndef __processorPlus_H__
#define __processorPlus_H__ 

#include "mpiTypesPlus.hpp"

namespace pFlow::Plus
{
	bool checkMPI(
		const char* funcName,
		int error, 
		bool forceAbort, 
		const char* fileName, 
		int lineNumebr);

	int  parReportAndExit(int errorCode);
}

#define ExitParCode(errorCode)\
	pFlow::Plus::parReportAndExit((errorCode))

#ifndef CheckMPI
#define CheckMPI(caller, fAbort)\
   pFlow::Plus::checkMPI(#caller, (caller), fAbort, __FILE__, __LINE__);
#endif
   
namespace pFlow::Plus
{



class processor
{
public:
	
	static inline 
	bool isSelfInitialized_ = false;

	static 
	void initMPI(int argc, char *argv[]);

	static
	void finalizeMPI();

	processor();
	
	~processor();

	static
	int myProcessorNo();
	
	static
	int nProcessors();
	
	static
	int masterNo();

	static
	bool isParallel();
	
	static
	bool isInitialized();
	
	static
	bool isFinalized();
	
	static
	bool isMaster();
	
	static
	Comm worldCommunicator();
	
	static
	int commSize();
	
	static
	int commRank();

	static
	void abort(int error);

}; //processor



}



#endif //__processor_H__

