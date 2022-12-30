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

#ifndef __scatteredCommunication_hpp__ 
#define __scatteredCommunication_hpp__

#include "span.hpp"
#include "procCommunication.hpp"
#include "mpiCommunication.hpp"
#include "procVector.hpp"

namespace pFlow::MPI
{

template<typename T>
void performSum(span<T>& dest, T* src, span<const int32>& map);


template<typename T>
class scatteredCommunication
:
	public procCommunication
{
protected:

	bool 							indexCreated_ = false;

	procVector<span<const int32>> 	dataMaps_ {true};

	procVector<DataType>			indexedMap_ {true};

	procVector<uniquePtr<T>>		buffers_ {true};

	procVector<size_t> 				buffersSize_ {0, true};

	bool createIndexTypes()
	{

		if(indexCreated_)
		{
			for(auto& map:indexedMap_)
			{
				bool res = CheckMPI(MPI_Type_free(&map), false );
				if(!res) return false;
			}	
		}

		procCommunication proc;

		for(size_t i = 0; i<dataMaps_.size(); i++)
		{
			
			if(auto [newT, success] =  
				proc.createIndexedDataType<T>(dataMaps_[i]); success)
			{
				indexedMap_[i] = newT;
			}
			else
			{
				return false;
			}
		}

		indexCreated_ = true;
		
		return true;
	}

	bool checkForBuffers()
	{
		for(size_t i = 0; i< buffers_.size(); i++)
		{
			if(dataMaps_[i].size() > buffersSize_[i])
			{
				size_t newSize = (1.1 * dataMaps_[i].size())+1;
				buffers_[i].reset(new T[newSize]);
				buffersSize_[i] = newSize;
			}
		}
		return true;
	}

public:

	scatteredCommunication()= default;

	scatteredCommunication(procVector<span<const int32>>& maps)
	{
		if(!changeDataMaps(maps))
		{
			processor::abort(0);
		}	
	}

	~scatteredCommunication()=default;

	scatteredCommunication(const scatteredCommunication&)=delete;

	scatteredCommunication& operator=(const scatteredCommunication&) = delete;

	bool changeDataMaps(procVector<span<const int32>>& maps)
	{
		

		if(processor::isMaster())
		{
			
			dataMaps_ = maps;
			if(!createIndexTypes())
			{
				fatalErrorInFunction<<"failed to create index types "<<endl;
				return false;
			}

			if(!checkForBuffers())
			{
				fatalErrorInFunction<<
				"failed to allocate buffer for scatteredCommunication"<<endl;
				return false;
			}
		}
		return true;
	}

	bool distribute(span<T>& sendBuff, span<T>& recvb)
	{
		
		procVector<Request> requests(true);
		procVector<Status> statuses(true);

		if(processor::isMaster())
		{
			bool res = true;
			for(size_t i = 0; i<indexedMap_.size(); i++)
			{
				res = res&&CheckMPI(
					MPI_Isend( 
						sendBuff.data(), 
						1, 
						indexedMap_[i], 
						i, 
						0, 
						processor::worldCommunicator(),
						&requests[i]), 
					false);
			}

			if(!res)return false;		
		}

		
		MPI_Status stat;	
		bool sucss = CheckMPI( 
			MPI_Recv(
				recvb.data(), 
				recvb.size()*sFactor<T>(), 
				Type<T>(), 
				0, 
				0, 
				processor::worldCommunicator(),
				&stat),
			false);


		if(processor::isMaster())
		{
			CheckMPI(
				MPI_Waitall(requests.size(), requests.data(), statuses.data()),
				false
				);
		}

		return sucss;
	}

	bool collectSum(span<T>& sendBuff, span<T>& recvb)
	{
		
		bool succs;
	
		succs = CheckMPI(
			MPI_Isend(
				sendBuff.data(),
				sendBuff.size()*sFactor<T>(),
				Type<T>(),
				0,
				0,
				processor::worldCommunicator(),
				&RequestNull),
			false);
		if(!succs) return false;

		procVector<Request> requests;

		if(processor::isMaster())
		{
			succs = true;
			for(size_t i=0; i<buffers_.size(); i++)
			{				
				succs = succs&& CheckMPI(
					MPI_Irecv(
						buffers_[i].get(),
						buffersSize_[i]*sFactor<T>(),
						Type<T>(),
						i,
						0,
						processor::worldCommunicator(),
						&requests[i]
						),
					true);
			}

			if(!succs) return false;
			Status stat;
			size_t numFinished = 0;
			int    procNo;
			while (numFinished < requests.size())
			{
				MPI_Waitany(requests.size(), requests.data(), &procNo, &stat);
				numFinished++;
				performSum(recvb, buffers_[procNo].get(), dataMaps_[procNo]);
			}

		}

		return true;
	}

};


template<typename T>
void performSum(span<T>& dest, T* src, span<const int32>& map)
{
	for(size_t i=0; i<map.size(); i++)
	{
		dest[map[i]] += src[i]; 
	}
}


} //pFlow::MPI

#endif //__scatteredCommunication_hpp__
