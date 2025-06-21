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

#ifndef __adaptiveGaussian_hpp__
#define __adaptiveGaussian_hpp__

#include "distribution.hpp"

namespace pFlow::coupling
{

class couplingMesh;

class adaptiveGaussian
:
    public distribution
{
    
   // radius of circule for cell neighbor list 
	Foam::scalar 				distLength_;

    Foam::label                 maxLayers_;

    const static inline Foam::scalar   a_ = 0.6142275;
    const static inline Foam::scalar   exponent_ = -0.6195039;

    inline 
    Foam::scalar stdDeviation(const Foam::scalar& dx_dp, const Foam::scalar& dcell) const
    {
        return dcell* a_ * Foam::pow(dx_dp, exponent_);
    }

public:

    /// Type info
    TypeInfoNV("adaptiveGaussian");

    /// Construc from dictionary 
    adaptiveGaussian(
        Foam::dictionary 		dict, 
        const couplingMesh& 	cMesh,
        const Plus::centerMassField& centerMass);

    /// Destructor
    ~adaptiveGaussian() = default;

    void updateWeights
    (
        const Plus::procCMField<Foam::label> & parCellIndex,
        const Plus::procCMField<real> & parDiameter
    );

}; // end adaptiveGaussian

} // end pFlow::coupling


#endif //__adaptiveGaussian_hpp__
