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

#ifndef __mpiCommunicationPlus_H__
#define __mpiCommunicationPlus_H__


#include "mpiTypesPlus.hpp"
#include "span.hpp"


namespace pFlow::Plus
{

extern DataType realx3Type__;

extern DataType int32x3Type__;


template<typename T> 
inline
auto constexpr Type()
{
	return MPI_BYTE;
}

template<typename T>
inline
auto constexpr sFactor()
{
	return sizeof(T);
}

template<char> 
inline
auto constexpr Type()
{
	return MPI_CHAR;
}
template<char>
inline
auto constexpr sFactor()
{
	return 1;
}

template<short>
auto constexpr Type()
{
	return MPI_SHORT;
}
template<short>
inline
auto constexpr sFactor()
{
	return 1;
}

template<unsigned short>
auto constexpr Type()
{
	return MPI_UNSIGNED_SHORT;
}
template<unsigned short>
auto constexpr sFactor()
{
	return 1;
}

template<int>
auto constexpr Type()
{
	return MPI_INT;
}
template<int>
auto constexpr sFactor()
{
	return 1;
}

template<>
auto constexpr Type<unsigned int>()
{
	return MPI_UNSIGNED;
}
template<>
auto constexpr sFactor<unsigned int>()
{
	return 1;
}

template<>
auto constexpr Type<long>()
{
	return MPI_LONG;
}
template<>
auto constexpr sFactor<long>()
{
	return 1;
}

template<>
auto constexpr Type<unsigned long>()
{
	return MPI_UNSIGNED_LONG;
}
template<>
auto constexpr sFactor<unsigned long>()
{
	return 1;
}


template<>
auto constexpr Type<float>()
{
	return MPI_FLOAT;
}
template<>
auto constexpr sFactor<float>()
{
	return 1;
}

template<>
auto constexpr Type<double>()
{
	return MPI_DOUBLE;
}
template<>
auto constexpr sFactor<double>()
{
	return 1;
}

template<>
inline auto Type<realx3>()
{
	return realx3Type__;
}
template<>
inline auto constexpr sFactor<realx3>()
{
	return 1;
}

template<>
inline auto Type<int32x3>()
{
	return int32x3Type__;
}
template<>
inline auto constexpr sFactor<int32x3>()
{
	return 1;
}

inline auto TypeCommit(DataType* type)
{
	return MPI_Type_commit(type);
}


template<typename T>
inline auto getCount(Status* status, int& count)
{
	int lCount;
	auto res = MPI_Get_count(status, Type<T>(), &lCount);
	count = lCount/sFactor<T>();
	return res;
}

template<typename T>
inline int convertIndex(const int& ind)
{
	return ind*sFactor<T>();
}

template<typename T> 
inline auto send(T* data, int count, int dest, int tag, Comm comm)
{
	return MPI_Send(data, sFactor<T>()*count, Type<T>(), dest, tag, comm);
}

template<typename T> 
inline auto send(T* data, int count, int dest, Comm comm)
{
	return MPI_Send(data, sFactor<T>()*count, Type<T>(), dest, 0, comm);
}


template<typename T>
inline auto recv(T* data, int count, int source, int tag, Comm comm, Status *status)
{
	return MPI_Recv(data, sFactor<T>()*count, Type<T>(), source, tag, comm, status);
}

template<typename T>
inline auto recv(T* data, int count, int source, Comm comm, Status *status)
{
	return MPI_Recv(data, sFactor<T>()*count, Type<T>(), source, 0, comm, status);
}

template<typename T>
inline auto scan(T sData, T& rData, Comm comm, Operation op = SumOp)
{
	return MPI_Scan(&sData, &rData, sFactor<T>()*1, Type<T>(), op , comm );
}

// gathering one scalar data to root processor 
template<typename T>
inline auto gather(T sendData, span<T>& recvData, int root, Comm comm)
{
	return MPI_Gather(
		&sendData, 
		sFactor<T>()*1, 
		Type<T>(), 
		recvData.data(),
		sFactor<T>()*1,
		Type<T>(),
		root,
		comm);
}

template<typename T>
inline auto allGather(T sendData, span<T>& recvData, Comm comm)
{
	return MPI_Allgather(
		&sendData,
		sFactor<T>()*1,
		Type<T>(),
		recvData.data(),
		sFactor<T>()*1,
		Type<T>(),
		comm);
}

template<typename T>
inline auto scatter(span<T>& sendData, T& recvData, int root, Comm comm)
{
	return MPI_Scatter(
		sendData.data(),
		sFactor<T>()*1,
		Type<T>(),
		&recvData,
		sFactor<T>()*1,
		Type<T>(),
		root,
		comm);
}

template<typename T>
inline auto Bcast(T& sendData, int root, Comm comm)
{
	return MPI_Bcast(
		&sendData, sFactor<T>()*1, Type<T>(), root, comm);

}

inline auto Wait(Request* request, Status* status)
{
	return MPI_Wait(request, status);
}

inline auto typeFree(DataType& type)
{
	return MPI_Type_free(&type);
}

// FILE operations 

inline auto fileOpen(Comm comm, const char* fileName, int amode, File& fh)
{
	return MPI_File_open(comm, fileName, amode, InfoNull, &fh);
}

inline auto fileClose(File& fh)
{
	return MPI_File_close(&fh);
} 

template<typename T>
inline auto fileWriteAt(File& fh, Offset off, T data, Status* status)
{
	return MPI_File_write_at(fh, off, &data, sFactor<T>()*1, Type<T>(), status);
}

template<typename T>
inline auto fileWriteAt(File& fh, Offset off, const span<T>& data, Status* status)
{
	return MPI_File_write_at(
		fh,
		off,
		data.data(),
		sFactor<T>()*data.size(),
		Type<T>(),
		status);
}

template<typename T>
inline auto fileWriteAtAll(File& fh, Offset off, T data, Status* status)
{
	return MPI_File_write_at_all(fh, off, &data, sFactor<T>()*1, Type<T>(), status);
}

template<typename T>
inline auto fileWriteAtAll(File& fh, Offset off, const span<T>& data, Status* status)
{
	return MPI_File_write_at_all(
		fh,
		off,
		data.data(),
		sFactor<T>()*data.size(),
		Type<T>(),
		status);
}

template<typename T>
inline auto fileReadAt(File& fh, Offset off, T& data, Status* status)
{
	return MPI_File_read_at(
		fh,
		off,
		&data,
		sFactor<T>()*1,
		Type<T>(),
		status);
}

template<typename T>
inline auto fileReadAt(File& fh, Offset off, span<T>& data, Status* status)
{
	return MPI_File_read_at(
		fh,
		off,
		data.data(),
		sFactor<T>()*data.size(),
		Type<T>(),
		status);
}

template<typename T>
inline auto fileReadAtAll(File& fh, Offset off, T& data, Status* status)
{
	return MPI_File_read_at_all(
		fh,
		off,
		&data,
		sFactor<T>()*1,
		Type<T>(),
		status);
}

template<typename T>
inline auto fileReadAtAll(File& fh, Offset off, span<T>& data, Status* status)
{
	return MPI_File_read_at_all(
		fh,
		off,
		data.data(),
		sFactor<T>()*data.size(),
		Type<T>(),
		status);
}




template<typename T>
inline auto fileIWriteAtAll(File& fh, Offset off, T data, Request* request)
{
	return MPI_File_iwrite_at_all(fh, off, &data, sFactor<T>()*1, Type<T>(), request);
}

template<typename T>
inline auto fileIWriteAtAll(File& fh, Offset off, const span<T>& data, Request* request)
{
	return MPI_File_iwrite_at_all(
		fh,
		off,
		data.data(),
		sFactor<T>()*data.size(),
		Type<T>(),
		request);
}




}


#endif  //__mpiCommunicationPlus_H__
