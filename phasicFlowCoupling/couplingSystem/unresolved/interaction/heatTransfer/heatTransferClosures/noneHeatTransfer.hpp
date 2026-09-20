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

#ifndef pFlow_noneHeatTransfer_hpp 
#define pFlow_noneHeatTransfer_hpp

#include "OFCompatibleHeader.hpp"
#include "typeInfo.hpp"

namespace pFlow
{
namespace coupling
{

class noneHeatTransfer
{
public:

    //- Type info
    
        TypeInfoNV("none");

    //- constructors

        explicit noneHeatTransfer(const Foam::dictionary& dict);

    //- public methods

        inline
        Foam::scalar dimlessNusselt(
            Foam::scalar Re, 
            Foam::scalar Pr, 
            Foam::scalar ep) const
        {
            return 0.0;
        }
        
        inline 
        Foam::scalar operator()(
            Foam::scalar Re, 
            Foam::scalar Pr, 
            Foam::scalar ep) const
        {
            return dimlessNusselt(Re, Pr, ep);
        }
        
}; 

} // coupling
} // pFlow

#endif // pFlow_noneHeatTransfer_hpp


