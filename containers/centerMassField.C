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

#include "centerMassField.hpp"


pFlow::MPI::centerMassField::centerMassField(size_t size, size_t capacity)
:
	eventSubscriber(),
	std::vector<realx3>()
{
	this->reserve(capacity);
	this->resize(size);
}

bool pFlow::MPI::centerMassField::checkForNewSize(size_t newSize)
{	

		if(newSize == this->size()) return true;
		
		eventMessage msg;
		if(newSize < this->capacity() )
		{
			// enough space is avaiable
			// only resize the container
			msg.add(eventMessage::SIZE_CHANGED);
			this->resize(newSize);
		}
		else if(newSize > this->capacity())
		{
			// resize to new size and let std::vector
			// decides about capacity 
			msg.add(eventMessage::SIZE_CHANGED);
			msg.add(eventMessage::CAP_CHANGED);
			this->resize(newSize);
		}

		if( !msg.isNull() )
		{
			return this->notify(msg);
		}

		return true;

}