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
#ifndef __procCommunicationPlus_hpp__
#define __procCommunicationPlus_hpp__

#include "processorPlus.hpp"
#include "procVectorPlus.hpp"
#include "mpiCommunicationPlus.hpp"

namespace pFlow::Plus
{


class procCommunication
:
	public processor
{
protected:

public:
	
	procCommunication();
	
	~procCommunication()=default;

	// - send a single val to all processors including itself
	template<typename T>
	std::pair<T,bool> distributeMasterToAll(const T& val)
	{
		
		T retVal = val;
		auto res = CheckMPI(Bcast(retVal, masterNo(),worldCommunicator()),false);
		
		return {retVal, res};
		
	}

	template<typename T>
	bool distributeMasterToAll(const T& val, T& recvVal)
	{
		recvVal = val;
		return CheckMPI(Bcast(recvVal, masterNo(),worldCommunicator()),false);
	}

	// send each values in vector (size is equal to number of processors) to each processor
	template<typename T>
	std::pair<T,bool> distributeMasterToAll(procVector<T>& vals)
	{

		T val;	
		auto vec = makeSpan(vals);
		auto res = CheckMPI(scatter(vec, val, masterNo(), worldCommunicator()),false);
		
		return {val, res};
	}

	template<typename T>
	std::pair<procVector<T>, bool> collectAllToAll(const T& val)
	{
		procVector<T> allVec;
		auto vec = makeSpan(allVec);
		auto res = CheckMPI(allGather(val, vec, worldCommunicator()), false);
		
		return {allVec, res};
	}

	template<typename T>
	bool collectAllToAll(const T& val, procVector<T>& allVec)
	{
		auto vec = makeSpan(allVec);
		return CheckMPI(allGather(val, vec, worldCommunicator()), false);
	}

	template<typename T>
	std::pair<procVector<T>,bool> collectAllToMaster(const T& val)
	{
		// only on master processor
		procVector<T> masterVec(true);
		
		auto masterSpan = makeSpan(masterVec);
		auto res = CheckMPI( 
			gather(val,masterSpan, masterNo(), worldCommunicator()), false);

		return {masterVec, res};

	}


	template<typename T> 
	std::pair<DataType,bool> createIndexedDataType(span<const int32> index)
	{
		DataType newType;
		auto res = CheckMPI( MPI_Type_create_indexed_block(
			index.size(), 
			sFactor<T>(), 
			index.data(), 
			Type<T>(), 
			&newType), false);
		
		if(res)
		{
			TypeCommit(&newType);
		}

		return {newType,res};	
	}

	template<typename T> 
	std::pair<DataType,bool> createIndexedDataType(const std::vector<int32>& index)
	{
		return createIndexedDataType<T>( span(index.data(), index.size()) );
	}



}; //procCommunication

} // pFlow::Plus

#endif //__procCommunicationPlus_hpp__
