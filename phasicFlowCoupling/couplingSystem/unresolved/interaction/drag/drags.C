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


#include "drag.hpp"
#include "self.hpp"
#include "Gaussian.hpp"
#include "Gaussian2.hpp"
#include "GaussianIntegral.hpp"
#include "adaptiveGaussian.hpp"
//#include "Gaussian29.hpp"

template class pFlow::coupling::drag<pFlow::coupling::self>;
template class pFlow::coupling::drag<pFlow::coupling::Gaussian>;
template class pFlow::coupling::drag<pFlow::coupling::Gaussian2>;
template class pFlow::coupling::drag<pFlow::coupling::GaussianIntegral>;
template class pFlow::coupling::drag<pFlow::coupling::adaptiveGaussian>;
//template class pFlow::coupling::drag<pFlow::coupling::Gaussian29>;
