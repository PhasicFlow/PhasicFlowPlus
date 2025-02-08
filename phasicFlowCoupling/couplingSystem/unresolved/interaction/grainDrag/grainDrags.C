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

#include "grainDrag.hpp"
#include "self.hpp"
#include "Gaussian.hpp"
#include "GaussianIntegral.hpp"

#include "DiFelice.hpp"
#include "Rong.hpp"
#include "ErgunWenYu.hpp"

#define makeGrainDrag(closure)                     \
													\
template class pFlow::coupling::grainDrag 			\
	< 												\
		pFlow::coupling::self, 						\
		closure, 									\
		true 										\
	>;  											\
template class pFlow::coupling::grainDrag  	 		\
	<												\
		pFlow::coupling::self, 						\
		closure, 									\
		false 										\
	>;  											\
													\
template class pFlow::coupling::grainDrag 			\
	<												\
		pFlow::coupling::Gaussian,					\
		closure,									\
		true 										\
	>; 												\
template class pFlow::coupling::grainDrag 			\
	< 												\
		pFlow::coupling::Gaussian, 					\
		closure,									\
		false										\
	>; 												\
template class pFlow::coupling::grainDrag 			\
	<												\
		pFlow::coupling::GaussianIntegral,			\
		closure,									\
		true										\
	>; 												\
template class pFlow::coupling::grainDrag 			\
	<												\
		pFlow::coupling::GaussianIntegral,			\
		closure,									\
		false										\
	>; 												

makeGrainDrag(pFlow::coupling::DiFelice);
makeGrainDrag(pFlow::coupling::ErgunWenYu);
makeGrainDrag(pFlow::coupling::Rong);
