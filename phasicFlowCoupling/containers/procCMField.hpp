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

#ifndef __procCMField_hpp__ 
#define __procCMField_hpp__

#include <vector>


#include "eventObserver.hpp"
#include "span.hpp"
#include "centerMassField.hpp"


namespace pFlow::MPI
{

template<typename T>
class procCMField
:
	public eventObserver,
	public std::vector<T>
{
protected:

	const centerMassField& centerMass_;

	bool 	preserveContent_ = false;

	word  fieldName_ = "procCMField_NoName";


	// restrict access to these methods 
	using std::vector<T>::reserve;

	using std::vector<T>::resize;

	using std::vector<T>::assign;

	using std::vector<T>::clear;

	using std::vector<T>::erase;

public:


	TypeInfoTemplate("procCMField",T);

	
	procCMField(const centerMassField& cm, bool subscribe = true)
	:
		eventObserver(cm, subscribe),
		std::vector<T>(cm.size(), cm.capacity()),
		centerMass_(cm)
	{}

	procCMField(
		word 	name, 
		const centerMassField& cm, 
		bool 	preserveContent = false, 
		bool 	subscribe = true)
	:
		procCMField(cm, subscribe)
	{
		fieldName_ = name;
		preserveContent_ = preserveContent;
	}

	procCMField(
		word name,
		const T& value,
		const centerMassField& cm,
		bool 	preserveContent = false,
		bool 	subscribe = true)
	:
		procCMField(name, cm, preserveContent, subscribe)
	{
		std::fill(this->begin(), this->end(), value);
	}


	procCMField(const procCMField&) = default;

	procCMField& operator = (const procCMField&) = default;

	procCMField& operator = (const std::vector<T>& src)
	{
		this->assign(src.begin(), src.end());
		return *this;
	}

	procCMField& operator =(const T& value)
	{
		std::fill(this->begin(), this->end(), value);
		return *this;
	}

	virtual ~procCMField()=default;

	
	bool update(const eventMessage& msg) override
	{
		// first check for capacity change 
		if(msg.isCapacityChanged())
		{
			if(!preserveContent_)
			{
				this->clear();
				this->reserve(centerMass_.capacity());
				this->resize(centerMass_.size());
			}else
			{
				this->reserve(centerMass_.capacity());
				this->resize(centerMass_.size());
			}

			return true;
		}

		if(msg.isSizeChanged())
		{
			this->resize(centerMass_.size());
			return true;
		}

		return true;
	}
	
};

template<typename T>
span<const T> makeSpan(const procCMField<T>& cmField)
{
	return span<const T>(cmField.data(), cmField.size());
}

template<typename T>
span<T> makeSpan(procCMField<T>& cmField)
{
	return span<T>(cmField.data(), cmField.size());
}


} // pFlow::MPI


#endif //__procCMField_hpp__
