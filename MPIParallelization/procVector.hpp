#ifndef __procVector_hpp__ 
#define __procVector_hpp__

// from std
#include <vector>

// from PhasicFlow
#include "iIstream.hpp"
#include "iOstream.hpp"
#include "processor.hpp"
#include "span.hpp"



namespace pFlow::MPI
{

template<typename T>
class procVector
:
	public std::vector<T>
{
public:

	using ProcVectorType = procVector<T>;

	using VectorType = std::vector<T>;

protected:

	int32 rank_ 	= 0;

	using VectorType::reserve;

	using VectorType::resize;

	using VectorType::assign;

	using VectorType::clear;

	using VectorType::erase;

public:


	procVector(bool onlyMaster = false)
	:
	VectorType(),
	rank_(pFlow::MPI::processor::myProcessorNo())
	{
		if( onlyMaster && !pFlow::MPI::processor::isMaster() ) return;

		this->reserve(pFlow::MPI::processor::nProcessors());
		this->resize(pFlow::MPI::processor::nProcessors());
	}

	procVector(const T& val, bool onlyMaster=false)
	:
	rank_(pFlow::MPI::processor::myProcessorNo())
	{
		if( onlyMaster && !pFlow::MPI::processor::isMaster() ) return;

		this->reserve(pFlow::MPI::processor::nProcessors());
		this->assign(pFlow::MPI::processor::nProcessors(),val);
	}

	procVector(const procVector&) = default;
	
	procVector(procVector&&) = default;

	procVector& operator=(const procVector&) = default;
	
	procVector& operator=(procVector&&) = default;

	procVector(const VectorType& src)
	:
		VectorType(src),
		rank_(pFlow::MPI::processor::myProcessorNo())
	{
		if(src.size()!= pFlow::MPI::processor::nProcessors())
		{
			fatalErrorInFunction<<
			"Size of std::vector and procVector does not match in copy"<<endl;
			processor::abort(0);
		}
	}

	procVector& operator=(const VectorType& src)
	{
		if(src.size() != this->size())
		{
			fatalErrorInFunction<<
			"Size of std::vector and procVector does not match in copy assignment"<<endl;
			processor::abort(0);
		}

		static_cast<VectorType&>(*this).operator=(src);
		return *this;
	}

	procVector& operator=(VectorType&& src)
	{
		if(src.size() != this->size())
		{
			fatalErrorInFunction<<
			"Size of std::vector and procVector does not match in move assignment"
			<<endl;
			processor::abort(0);
		}

		static_cast<VectorType&>(*this).operator=(std::move(src));
		return *this;
	}

	procVector(VectorType&& src)
	:
		VectorType(std::move(src)),
		rank_(pFlow::MPI::processor::myProcessorNo())
	{
		if(this->size()!= 
			static_cast<size_t>(pFlow::MPI::processor::nProcessors()))
		{
			fatalErrorInFunction<<
			"Size of std::vector and procVector does not match in move"<<endl;
			processor::abort(0);
		}
	}

	~procVector()=default;

	auto& thisValue()
	{
		return VectorType::operator[](rank_);
	}

	const auto& thisValue()const
	{
		return VectorType::operator[](rank_);
	}

	auto commSize()const
	{
		return VectorType::size();
	}

	auto commRank()const
	{
		return rank_;
	}

	auto worldCommunicator()const
	{
		return pFlow::MPI::processor::worldCommunicator();
	}

	auto worldCommunicator()
	{
		return pFlow::MPI::processor::worldCommunicator();
	}

	bool write(iOstream& os)const
	{
		// start
		os << token::BEGIN_LIST;
		for(const auto& vi: *this)
		{
			os << vi;
			if(vi != this->end())
				os<< token::NL;
		}

		os << token::END_LIST;

	    os.check(FUNCTION_NAME);
	    return true;  
	}

};

template<typename T>
span<T> makeSpan(procVector<T>& vec)
{
	return span(vec.data(), vec.size());
}

template<typename T>
span<const T> makeSpan(const procVector<T>& vec)
{
	return span(vec.data(), vec.size());
}

/*template<typename T>
inline iIstream& operator >> (iIstream & is, procVector<T> & ivec )
{
	if( !ivec.read(is) )
	{
		ioErrorInFile (is.name(), is.lineNumber());
		fatalExit;
	}
	return is;
}*/

template<typename T> 
inline iOstream& operator << (iOstream& os, const procVector<T>& ovec )
{
	
	if( !ovec.write(os) )
	{
		ioErrorInFile(os.name(), os.lineNumber());
		fatalExit;
	}

	return os; 
}

}

#endif
