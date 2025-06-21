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

#ifndef __fluidVelocity_hpp__
#define __fluidVelocity_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class couplingMesh;

class fluidVelocity
{
    const Foam::volVectorField&     Uf_;

    const bool                      interpolate_ = false;

    Plus::realx3ProcCMField         Up_;

public:

    fluidVelocity(
        const word& type, 
        const Foam::volVectorField& U,
        const couplingMesh& cMesh);

    void interpolate(const couplingMesh& cMesh);

    inline
    Foam::vector uFluid(Foam::scalar celli, size_t parIdx)const
    {
        /*if(interpolate_)
        {
            const auto& up = Up_[parIdx]; 
            return Foam::vector{up.x(), up.y(), up.z()};
        }
        else
        {
            return Uf_[celli];
        }*/

        return Uf_[celli];
    }

};

}


#endif //__fluidVelocity_hpp__
