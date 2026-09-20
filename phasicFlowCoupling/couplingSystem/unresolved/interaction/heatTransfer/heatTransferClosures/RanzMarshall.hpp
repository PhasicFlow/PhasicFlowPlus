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

#ifndef pFlow_RanzMarshall_hpp 
#define pFlow_RanzMarshall_hpp

#include "OFCompatibleHeader.hpp"
#include "typeInfo.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class RanzMarshall
 * @brief Implementation of the Ranz-Marshall heat transfer closure model.
 *
 * @details
 * Calculates the dimensionless Nusselt number (Nu) for spherical 
 * particles based on the local Reynolds (Re) and Prandtl (Pr) numbers.
 * The Ranz-Marshall correlation is highly suitable for isolated 
 * particles or highly dilute flows.
 * 
 * Nu = 2.0 + 0.6 * Re^(1/2) * Pr^(1/3)
 */
class RanzMarshall
{
public:

    //- Type info

        TypeInfoNV("RanzMarshall");

    //- constructors

        /**
         * @brief Constructor initializing the model from the simulation dict.
         * @param dict Sub-dictionary containing specific model coefficients.
         */
        explicit RanzMarshall(const Foam::dictionary& dict);

        ~RanzMarshall() = default;

    //- public methods

        /**
         * @brief Calculates the dimensionless Nusselt number.
         * @param Re Particle Reynolds number [-].
         * @param Pr Fluid Prandtl number [-].
         * @param ep Fluid volume fraction (porosity) [-]. Passed for 
         *           interface consistency but not utilized here.
         * @return The computed Nusselt number [-].
         */
        inline
        Foam::scalar dimlessNusselt(
            Foam::scalar Re, 
            Foam::scalar Pr, 
            Foam::scalar ep) const
        {
            return 2.0 + 0.6 * Foam::sqrt(Re) * Foam::pow(Pr, 1.0/3.0);
        }
        
        /**
         * @brief Functor operator for direct function-like calls.
         */
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

#endif // pFlow_RanzMarshall_hpp


