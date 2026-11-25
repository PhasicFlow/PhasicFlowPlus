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

#include "sphereMassTransfer.hpp"

// Nusselt closures of heatTransferModel, used here for the Sherwood number
#include "RanzMarshall.hpp"
#include "Gunn.hpp"
#include "noneHeatTransfer.hpp"

namespace pFlow
{

namespace coupling
{

    template class sphereMassTransfer<RanzMarshall>;

    template class sphereMassTransfer<Gunn>;

    template class sphereMassTransfer<noneHeatTransfer>;

} // coupling

} // pFlow
