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
 * @class Beetstra
 * @brief Drag closure model for intermediate Reynolds number flows in granular systems
 * 
 * This class implements the Beetstra drag correlation for computing drag forces
 * in particle-fluid systems. The model is applicable for intermediate Reynolds number
 * flows past mono-disperse arrays of spheres.
 * 
 * @details
 * The dimensionless drag force is calculated using the following formula:
 * 
 * \f[
 * f = \frac{10(1-\varepsilon)}{\varepsilon^2} + \varepsilon^2(1 + 1.5\sqrt{1-\varepsilon}) + 
 * \frac{0.413}{24\varepsilon^2} \cdot \frac{(\frac{1}{\varepsilon} + 3\varepsilon(1-\varepsilon) + 8.4 \text{Re}^{-0.343})}{1 + 10^{3(1-\varepsilon)} \text{Re}^{-0.5 + 2(1-\varepsilon)}} \cdot \text{Re}
 * \f]
 * 
 * where:
 * - \f$ \varepsilon \f$ is the void fraction (gas/liquid volume fraction)
 * - \f$ \text{Re} \f$ is the relative Reynolds number
 * - A residual Reynolds number is enforced to prevent numerical issues at very low Re values
 * 
 * @reference
 * Beetstra, R., van der Hoef, M. A., & Kuipers, J. A. M. (2007).
 * "Drag force of intermediate Reynolds number flow past mono- and bidisperse
 * arrays of spheres." AIChE Journal, 53(3), 489-501.
 * 
 * @note A residual Reynolds number is used internally to prevent numerical instabilities
 * at very low Re values.
 */

#ifndef __Beetstra_hpp__ 
#define __Beetstra_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

#include "typeInfo.hpp"




namespace pFlow::coupling
{
class Beetstra
{

    Foam::scalar 		residualRe_;
     
    

public:

    // type info
    TypeInfoNV("Beetstra");

    Beetstra(const Foam::dictionary&	dict);

    ~Beetstra() = default;


    inline
    Foam::scalar dimlessDrag(Foam::scalar Re, Foam::scalar ep)const
    {
        const Foam::scalar ep_p = 1-ep;
        const Foam::scalar ep2 = ep*ep;
        Re = Foam::max(Re, residualRe_);

        return 10*ep_p/ep2 +
            ep2*(1+1.5*Foam::sqrt(ep_p)) +
            (0.413/24.0/ep2)*
                (
                    (1/ep + 3*ep*ep_p + 8.4*Foam::pow(Re,-0.343)) /
                    (1+ Foam::pow(10.0,3.0*ep_p)*Foam::pow(Re, -0.5 + 2*ep_p))
                )*Re;
    }

    inline 
    Foam::scalar operator()(Foam::scalar Re, Foam::scalar ep)const
    {
        return dimlessDrag(Re, ep);
    }
        
    
}; 

} // pFlow::coupling


#endif
