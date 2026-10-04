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

#ifndef __Cello_hpp__ 
#define __Cello_hpp__

#include "OFCompatibleHeader.hpp"

#include "typeInfo.hpp"

namespace pFlow::coupling
{

class Cello
{

    Foam::scalar residualRe_;

public:

    TypeInfoNV("Cello");

    Cello(const Foam::dictionary& dict);

    ~Cello() = default;

    inline
    Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep) const
    {
        const Foam::scalar Rec = Foam::max(Re, residualRe_);
        const Foam::scalar epc = Foam::max(ep, Foam::SMALL);
        const Foam::scalar ep2 = epc*epc;
        const Foam::scalar ep3 = ep2*epc;
        const Foam::scalar ep4 = ep2*ep2;
       

        const Foam::scalar K0 = (1.0 - epc)/(1.0 + 3.0*epc);
        const Foam::scalar K1 = (1.0 + 128.0*K0 + 715.0*Foam::sqr(K0))
            /(ep2*(1.0 + 49.5*K0));

        const Foam::scalar Re2 = Rec*Rec;

        const Foam::scalar K2 =
            (1.0 + 0.130*Rec + 6.66e-4*Re2)
            /(1.0 + 3.42e-2*Rec + 6.92e-6*Re2);

        const Foam::scalar K3n = -410.0*epc + 9.20e7*Rec*Foam::pow(K0,20)
            + 1900.0*ep2 - 6.60e-2*Rec;

        const Foam::scalar K3d = 6600.0*epc + 4.92e-4*Re - 4.30e4*ep2
            - 1.31e-4*Re2 + 7.38e4*ep3;

        const Foam::scalar K3 =
            (2.0*Re2/(1.0 + Rec))
            *(K3n/Foam::max(K3d, Foam::SMALL));

        return K1 + K2*ep4 + K3*(1.0 - ep4);
    }

    inline
    Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep) const
    {
        return dimlessDrag(Re, ep);
    }

};

} // pFlow::coupling

#endif

