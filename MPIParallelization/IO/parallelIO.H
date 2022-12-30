#ifndef __parallelIO_H__ 
#define __parallelIO_H__ 

#include "processors.H"
#include "Vector.H"
#include "fileSystem.H"
#include "span.H"

namespace dFlow::MPI
{


class parallelIO
:
	public processors
{
public:
	
	enum DataFormat: FileSize
	{
		Binary = 1,
		ASCII  = 2
	};

protected:
	
	Request 		ioRequest_ = RequestNull;

	File 			fh_ 		= FileNull;

	const bool 		blocking_ 	= true;

	inline 
	FileSize chunkSizeOffeset(int procNo)
	{
		// 2 is added for number of processors and data format flag
		return (2+procNo)*sizeof(FileSize);
	}

public:

	parallelIO(bool blocking);

	~parallelIO();

	bool waitForStreaming();
	
	bool writeToFile(const fileSystem& filePath, span<char> data, bool binary);

	// parallel reading of meta data 
	bool readMeta(const fileSystem& filePath, bool& binary, Vector<FileSize>& chunkSizes);

	// non-parallel reading of meta (only on master)
	bool readMetaMaster(const fileSystem& filePath, bool& binary, Vector<FileSize>& chunkSizes);

	bool readData(
		const fileSystem& filePath,
	 	bool binary, 
	 	const Vector<FileSize>& chunkSizes,
	 	span<char>& data);

	bool readDataMaster(
		const fileSystem& filePath,
	 	bool binary, 
	 	const Vector<FileSize>& chunkSizes,
	 	span<char>& data);

	
};

}

#endif //__parallelIO_H__ 
