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
 * @class ErgunWenYu
 * @brief Ergun-Wen-Yu drag correlation for fluid-particle
 *        interactions in unresolved coupling simulations.
 *
 * This class implements the Ergun-Wen-Yu drag model which combines
 * the Ergun equation for dense flow regimes (ε < 0.8) with the Wen-Yu
 * correlation for dilute flow regimes (ε ≥ 0.8).
 *
 * @section Formulation
 *
 * The dimensionless drag force is calculated as follows:
 *
 * For dilute flows (ε ≥ 0.8):
 * @f[
 * f = \frac{C_d}{24} Re \varepsilon^{-3.65}
 * @f]
 * where
 * @f[
 * C_d = \begin{cases}
 * 24(1 + 0.15 Re^{0.687})/Re & \text{if } Re \leq 1000 \\
 * 0.44 & \text{if } Re > 1000
 * \end{cases}
 * @f]
 *
 * For dense flows (ε < 0.8):
 * @f[
 * f = \frac{150(1-\varepsilon) + 1.75 Re}{18 \varepsilon^2}
 * @f]
 *
 * where Re is the particle Reynolds number and ε is the
 * fluid volume fraction.
 *
 * @references
 * - Gidaspow, D. (1994) Multiphase Flow and Fluidization:
 *   Continuum and Kinetic Theory Description, Academic Press,
 *   San Diego, CA.
 * - Ergun, S. (1952) Fluid flow through packed columns.
 *   Chemical Engineering Progress, 48, 89–94.
 * - Wen, C.Y. and Yu, Y.H. (1966) Mechanics of fluidization.
 *   Chemical Engineering Progress Symposium Series, 62, 100–111.
 */


#ifndef __ErgunWenYu_hpp__ 
#define __ErgunWenYu_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

#include "typeInfo.hpp"


namespace pFlow::coupling
{
class ErgunWenYu
{

	Foam::scalar 		residualRe_;
	 
	

public:

	// type info
	TypeInfoNV("ErgunWenYu");

	ErgunWenYu(const Foam::dictionary&	dict);

	~ErgunWenYu() = default;


	inline
	Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)const
	{
		if( ep >= 0.8 )
		{
			Foam::scalar Cd;
			if(Re <= 1000.0 )
				Cd = 24 * ( 1+0.15*Foam::pow(Re,0.687) ) / Re; 
			else
				Cd = 0.44;
			return Cd/24 * Re * Foam::pow(ep, -3.65 ); 
		}else
		{
			return 	(150.0*(1.0-ep)+ 1.75*Re )/(18.0*ep*ep);
		}

	}

	inline 
	Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
	{
		return dimlessDrag(Re, ep);
	}
		
	
}; 

} // pFlow::coupling


#endif
