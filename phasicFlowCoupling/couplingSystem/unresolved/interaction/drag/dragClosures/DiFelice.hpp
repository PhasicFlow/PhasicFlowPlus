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


/**
 * @class DiFelice
 * @brief Di Felice drag coefficient correlation for fluid-particle interaction systems.
 * 
 * This class implements the drag coefficient correlation proposed by Di Felice (1994)
 * for fluid-particle interaction systems. The correlation accounts for the effect of
 * void fraction (voidage) on the drag force through an exponential function.
 * 
 * @details
 * The dimensionless drag force is calculated using the Di Felice correlation:
 * 
 * \f[
 * \f = \frac{C_d}{24} \cdot Re \cdot \varepsilon^{-\xi}
 * \f]
 * 
 * where:
 * - \f$C_d = \left(0.63 + \frac{4.8}{\sqrt{Re}}\right)^2\f$ is the drag coefficient
 * - \f$Re\f$ is the Reynolds number
 * - \f$\varepsilon\f$ is the void fraction (volume fraction of fluid)
 * - \f$\xi = 3.7 - 0.65 \exp\left(-0.5(1.5 - \log_{10}Re)^2\right)\f$ is the voidage exponent
 * 
 * @reference
 * Di Felice, R. (1994) The voidage function for fluid–particle interaction systems.
 * International Journal of Multiphase Flow, 20, 153–159.
 * 
 * @see dimlessDrag()
 */
#ifndef __DiFelice_hpp__ 
#define __DiFelice_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"


#include "typeInfo.hpp"

namespace pFlow::coupling
{

class DiFelice
{

	Foam::scalar 		residualRe_;

public:

	// type info
	TypeInfoNV("DiFelice");

	DiFelice(const Foam::dictionary& dict);

	inline
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep) const
	{

		auto Rec = Foam::max(Re,residualRe_);
		Foam::scalar xi = 3.7 - 0.65*Foam::exp(-0.5*Foam::pow(1.5-Foam::log10(Rec),2));
		Foam::scalar Cd = Foam::pow(0.63+4.8/Foam::sqrt(Rec),2);

		return Cd/24 * Re * Foam::pow(ep, -xi ); 
		
	}
	
	inline 
	Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
	{
		return dimlessDrag(Re, ep);
	}
	
}; 

} // pFlow::coupling


#endif
