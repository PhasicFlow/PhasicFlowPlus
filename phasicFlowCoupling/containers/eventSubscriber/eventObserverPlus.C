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

#include "eventObserverPlus.hpp"
#include "eventSubscriberPlus.hpp" 

pFlow::Plus::eventObserver::eventObserver():
	subscriber_(nullptr),
	subscribed_(false)
{}

pFlow::Plus::eventObserver::eventObserver
(
	const eventSubscriber& subscriber,
	bool subscribe
)
:
	subscriber_(&subscriber),
	subscribed_(false)
{
	if(subscribe && subscriber_)
	{
		subscribed_ = subscriber_->subscribe(this);
	}
}

pFlow::Plus::eventObserver::~eventObserver()
{
	if(subscribed_ && subscriber_)
		subscriber_->unsubscribe(this);
}

bool pFlow::Plus::eventObserver::subscribe(const eventSubscriber& subscriber)
{
	subscriber_ = &subscriber;
	subscribed_ = subscriber_->subscribe(this);
	return subscribed_;
}

