#include <cstdio>

#include "parallelIO.H"
#include "mpiCommunication.H"
#include "error.H"
#include "streams.H"


 dFlow::MPI::parallelIO::parallelIO(bool blocking)
 :
 	processors(),
 	ioRequest_(RequestNull),
 	blocking_(blocking)
{

}

dFlow::MPI::parallelIO::~parallelIO()
{
	waitForStreaming();
}


bool dFlow::MPI::parallelIO::waitForStreaming()
{
	if(blocking_) return true;

	bool sB = CheckMPI(Wait(&ioRequest_, StatusIgnore), false);
	if ( sB )
	{
		if(fh_!= FileNull)
			return CheckMPI(fileClose(fh_), false);

		return true;
	}
	else
	{
		return false;
	}
}


bool dFlow::MPI::parallelIO::writeToFile(
	const fileSystem& filePath,
	span<char> data,
	bool binary)
{
	// collect information from other processes
	FileSize numProc = this->commSize();
	FileSize thisSize = data.size();
	FileSize offset;


	// scan the offset 
	CheckMPI(
		scan(thisSize, offset, worldCommunicator()), true);

	
	word wFile = filePath.wordPath();

	CheckMPI( 
		fileOpen(
			worldCommunicator(),
			wFile.c_str(),
			ModeCreate + ModeWriteOnly,
			fh_),
		true);

	// first location, number of processors
	if(isMaster())
	{
		CheckMPI(
			fileWriteAt(
				fh_, 
				0,
				(binary?DataFormat::Binary:DataFormat::ASCII),
				StatusIgnore),
			true);

		CheckMPI(fileWriteAt(fh_, sizeof(FileSize), numProc, StatusIgnore), true);
	}

	// data chunk size 
	FileSize sizeOffset = chunkSizeOffeset(this->commRank());
	CheckMPI(fileWriteAtAll(fh_, sizeOffset, thisSize, StatusIgnore), true);

	// the data chunk
	offset -= thisSize;
	FileSize chunkOffset = chunkSizeOffeset(this->commSize()) + offset;
	if(blocking_)
	{
		CheckMPI(fileWriteAtAll(fh_, chunkOffset, data, StatusIgnore), true);	
		CheckMPI( fileClose(fh_), true);
	}
	else
	{
		CheckMPI(fileIWriteAtAll(fh_, chunkOffset, data, &ioRequest_), true);	
	}
	

	

	return true;
}

bool dFlow::MPI::parallelIO::readMeta(
	const fileSystem& filePath,
	bool& binary,
	Vector<FileSize>& chunkSizes)
{
	
	word wFile = filePath.wordPath();
	
	CheckMPI( 
		fileOpen(
			worldCommunicator(),
			wFile.c_str(),
			ModeReadOnly,
			fh_),
		true);

	FileSize format, numProc;

	CheckMPI(
		fileReadAtAll(
			fh_,
			0,
			format,
			StatusIgnore),
		true);

	if(format == DataFormat::Binary)
	{
		binary = true;
	}
	else if(format == DataFormat::ASCII)
	{
		binary = false;
	}
	else
	{
		fatalErrorInFunction<<
		"Unrecognized DataFormat from file "<< filePath<<endl;
		return false;
	}

	CheckMPI(
		fileReadAtAll(
			fh_,
			sizeof(FileSize),
			numProc,
			StatusIgnore),
		true);

	if( numProc == 0 || numProc > MaxNoProcessors)
	{
		fatalErrorInFunction<<
		"Invalid number of processors read from file "<<filePath<<endl;
		return false;
	}

	chunkSizes.resize(numProc);

	auto spn = makeSpan(chunkSizes);

	CheckMPI(
		fileReadAtAll(
			fh_,
			2*sizeof(FileSize),
			spn,
			StatusIgnore),
		true);


	return true;
}

bool dFlow::MPI::parallelIO::readMetaMaster(
	const fileSystem& filePath,
	bool& binary,
	Vector<FileSize>& chunkSizes)
{

	// only on master
	if(!isMaster()) return true;

	word wFile = filePath.wordPath();
	auto fh = std::fopen(wFile.c_str(), "rb");

	if(!fh)
	{
		fatalErrorInFunction<<
		"Error in Opening file "<< filePath<<endl;
		return false;
	}

	FileSize format, numProc;

	auto res = std::fread(
			&format,
			sizeof(FileSize),
			1,
			fh);

	if(res != sizeof(FileSize))
	{
		fatalErrorInFunction<<
		"Error in reading file "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}
		
	if(format == DataFormat::Binary)
	{
		binary = true;
	}
	else if(format == DataFormat::ASCII)
	{
		binary = false;
	}
	else
	{
		fatalErrorInFunction<<
		"Unrecognized DataFormat from file "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}

	res = std::fread(&numProc, sizeof(FileSize), 1, fh);

	if(res != sizeof(FileSize))
	{
		fatalErrorInFunction<<
		"Error in reading file "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}

	if( numProc == 0 || numProc > MaxNoProcessors)
	{
		fatalErrorInFunction<<
		"Invalid number of processors read from file "<<filePath<<endl;
		std::fclose(fh);
		return false;
	}

	chunkSizes.resize(numProc);

	res = std::fread(chunkSizes.data(), sizeof(FileSize), numProc, fh);

	if(res!= numProc*sizeof(FileSize))
	{
		fatalErrorInFunction<<
		"Error in reading chunkSizes from file "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}
	std::fclose(fh);
	return true;
}


bool dFlow::MPI::parallelIO::readData(
	const fileSystem& filePath,
 	bool binary, 
 	const Vector<FileSize>& chunkSizes,
 	span<char>& data)
{

	if(chunkSizes.size() != this->commSize())
	{
		fatalErrorInFunction<< 
		"number of procssors ("<< yellowText(chunkSizes.size())<<
		") is not equal to number of chunks ("<<
		yellowText(this->commSize())<< ") in file "<< filePath<<endl;
		return false;
	}

	
	word wFile = filePath.wordPath();

	CheckMPI( 
		fileOpen(
			worldCommunicator(),
			wFile.c_str(),
			ModeReadOnly,
			fh_),
		true);

	auto procNo = this->commRank();
	FileSize toRecv = chunkSizes[procNo];

	// start of data chunks 
	FileSize offset = chunkSizeOffeset(this->commSize());

	for(auto i=0; i<procNo; i++)
	{
		offset += chunkSizes[i];
	}

	if( data.size() != toRecv )
	{
		fatalErrorInFunction<<
		"The span data size ("<<yellowText(data.size())<<
		") is not equalt to the requested chunk size ("<<
		yellowText(toRecv)<<").\n";
		return false;
	}

	Status status;

	CheckMPI( 
		fileReadAtAll(
			fh_,
			offset,
			data,
			&status),
		true);

	int recvCount;
	
	CheckMPI(
		getCount<char>(&status, recvCount), false);

	if( static_cast<FileSize>(recvCount) != toRecv)
	{
		fatalErrorInFunction;
		return false;
	}

	return true;

}

bool dFlow::MPI::parallelIO::readDataMaster(
		const fileSystem& filePath,
	 	bool binary, 
	 	const Vector<FileSize>& chunkSizes,
	 	span<char>& data)
{

	if(!isMaster()) return true;

	// sum of all chuncks
	FileSize toRecv = sum(chunkSizes);

	if( data.size() != toRecv )
	{
		fatalErrorInFunction<<
		"The span data size ("<<yellowText(data.size())<<
		") is not equalt to the requested chunk size ("<<
		yellowText(toRecv)<<").\n";
		return false;
	}

	word wFile = filePath.wordPath();
	auto fh = std::fopen(wFile.c_str(), "rb");

	if(!fh)
	{
		fatalErrorInFunction<<
		"Error in Opening file "<< filePath<<endl;
		return false;
	}

	// start of data chunks 
	FileSize offset = chunkSizeOffeset(chunkSizes.size());
	
	if(auto res = std::fseek(fh, offset, SEEK_SET); res!= 0 )
	{
		fatalErrorInFunction<<
		"Error in file seek "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}

	if(auto res = std::fread(
			data.data(),
			sizeof(char),
			data.size(),
			fh);
		res!= data.size() )
	{
		fatalErrorInFunction<<
		"Error in reading file "<< filePath<<endl;
		std::fclose(fh);
		return false;
	}

	std::fclose(fh);
	return true;
}