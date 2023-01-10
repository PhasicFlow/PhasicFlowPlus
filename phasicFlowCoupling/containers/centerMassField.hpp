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

#ifndef __centerMassField_hpp__ 
#define __centerMassField_hpp__

#include <vector>

#include "eventSubscriber.hpp"
#include "span.hpp"
#include "types.hpp"

namespace pFlow::MPI
{

class centerMassField
:
	public eventSubscriber,
	public std::vector<realx3>
{
protected:
	
	using std::vector<realx3>::reserve;

	using std::vector<realx3>::resize;

	using std::vector<realx3>::assign;

	using std::vector<realx3>::clear;

	using std::vector<realx3>::erase;

public:

	TypeInfo("centerMassField");

	centerMassField() = default;

	centerMassField(size_t size, size_t capacity);

	virtual ~centerMassField()=default;

	bool checkForNewSize(size_t newSize);
	
};

inline 
span<realx3> makeSpan(centerMassField& cmField)
{
	return span( cmField.data(), cmField.size());
}

inline 
span<const realx3> makeSpan(const centerMassField& cmField)
{
	return span<const realx3>( cmField.data(), cmField.size());
}


} // pFlow::MPI


#endif //__centerMassField_hpp__
