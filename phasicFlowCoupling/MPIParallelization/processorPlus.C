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
#include "Pstream.H"

// from PhasicFlow
#include "error.hpp"
#include "processorPlus.hpp"
#include "streams.hpp"
#include "mpiCommunicationPlus.hpp"


static int numVarsInitialized__ = 0;

pFlow::Plus::DataType pFlow::Plus::realx3Type__;

pFlow::Plus::DataType pFlow::Plus::int32x3Type__;


void pFlow::Plus::processor::initMPI(int argc, char *argv[])
{
	if(!processor::isInitialized())
	{
		CheckMPI(MPI_Init(&argc, &argv), true);
		isSelfInitialized_ = true;
	}
}

void pFlow::Plus::processor::finalizeMPI()
{
	if(isSelfInitialized_ && !isFinalized())
	{
		CheckMPI(MPI_Finalize(), true);
	}
}

pFlow::Plus::processor::processor()
{
	if(isParallel() && !isInitialized())
	{
		fatalErrorInFunction<<
		"MPI communication is not initialized yet!"<<endl;
		processor::abort(ErrOp);
	}

	if( numVarsInitialized__ == 0 ) 
	{
		MPI_Type_contiguous(3, Type<real>(), &realx3Type__);
		MPI_Type_commit(&realx3Type__);

		MPI_Type_contiguous(3, Type<int32>(), &int32x3Type__);
		MPI_Type_commit(&int32x3Type__);
	}

	numVarsInitialized__++;	
	
}

pFlow::Plus::processor::~processor()
{
	/*if(numVarsInitialized__==1)
	{
		MPI_Type_free(&realx3Type__);
		MPI_Type_free(&int32x3Type__);
		numVarsInitialized__ = 0;	
	}
	numVarsInitialized__ --;*/
}

int pFlow::Plus::processor::myProcessorNo()
{
	return Foam::UPstream::myProcNo();
}

int pFlow::Plus::processor::nProcessors()
{
	return Foam::UPstream::nProcs();
}

int pFlow::Plus::processor::masterNo()
{
	return 0;
}

bool pFlow::Plus::processor::isParallel()
{
	return Foam::UPstream::parRun();
}

bool pFlow::Plus::processor::isInitialized()
{
	int res;
	MPI_Initialized(&res);
	return res;
}

bool pFlow::Plus::processor::isFinalized()
{
	int res;
	MPI_Finalized(&res);
	return res;
}

bool pFlow::Plus::processor::isMaster()
{
	return Foam::UPstream::master();
}

pFlow::Plus::Comm pFlow::Plus::processor::worldCommunicator()
{
	return CommWorld;
}

int pFlow::Plus::processor::commSize()
{
	return Foam::UPstream::nProcs();
}

int pFlow::Plus::processor::commRank()
{
	return Foam::UPstream::myProcNo();
}

void pFlow::Plus::processor::abort(int error)
{
	MPI_Abort(processor::worldCommunicator(), error);	
}

bool pFlow::Plus::checkMPI(const char* funcName, int error, bool forceAbort, const char* fileName, int lineNumebr)
{
	if(error == MPI_SUCCESS) return true;
	output<< "Error occured in function "<< funcName <<
	" located in file "<< fileName<< " at line "<< lineNumebr<<endl;
	if(!forceAbort) return false;
	ExitParCode(error);
	return false;
}

int  pFlow::Plus::parReportAndExit(int errorCode)
{
	errReport<<"\n>>> dFlow is exiting . . ." << endl;
	processor::abort(errorCode);
	return errorCode;		
}


/*bool checkMPI(const char* funcName, int error, bool forceAbort, const char* fileName, int lineNumebr)
{
	if(error == MPI_SUCCESS) return true;
	pFlow::output<< "error occured in function "<< funcName <<
	" located in file "<< fileName<< " at line "<< lineNumebr<<pFlow::endl;
	if(!forceAbort) return false;
	ExitCode(error);
	return false;
}*/



