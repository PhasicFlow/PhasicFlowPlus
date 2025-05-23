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


#include "eventSubscriberPlus.hpp"
#include <algorithm>

pFlow::Plus::eventSubscriber::~eventSubscriber()
{
	for(auto& observer:observerList_)
	{
		observer->invalidateSubscriber();
	}
}

bool pFlow::Plus::eventSubscriber::subscribe
(
	eventObserver* observer
)const
{
	if(observer)
	{
		observerList_.push_back(observer);
		return true;
	}
	else
	{
		return false;
	}
}

bool pFlow::Plus::eventSubscriber::unsubscribe
(
	eventObserver* observer
)const
{
	if(observer)
	{
		std::remove(observerList_.begin(), observerList_.end(), observer);
	}
	return true;
}

bool pFlow::Plus::eventSubscriber::notify
(
	const eventMessage &msg
)
{
	for ( auto& observer:observerList_ )
	{
		if(observer)
			if( !observer->update(msg) ) return false;
	}

	return true;
}

/*bool pFlow::Plus::eventSubscriber::notify
(
	const eventMessage& msg,
	const List<eventObserver*>& exclutionList
)
{
	std::set<eventObserver*> sortedExcList(exclutionList.begin(),exclutionList.end());

	for(auto& observer:observerList_)
	{
		if( observer && sortedExcList.count(observer) == 0 )
		{
			if(!observer->update(msg)) return false;
		}
	}

	return true;
}*/