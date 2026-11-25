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

#ifndef pFlow_Gunn_hpp
#define pFlow_Gunn_hpp

#include "OFCompatibleHeader.hpp"
#include "typeInfo.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class Gunn
 * @brief Implementation of the Gunn heat transfer closure model.
 *
 * @details
 * Calculates the dimensionless Nusselt number (Nu) for spherical
 * particles based on the local Reynolds (Re) and Prandtl (Pr) numbers
 * and the fluid volume fraction (ep). In contrast to the Ranz-Marshall
 * correlation, which is derived for an isolated particle, the Gunn
 * correlation accounts for the presence of neighboring particles through
 * the porosity and is therefore suitable for fixed and fluidized beds.
 *
 * Nu = (7 - 10*ep + 5*ep^2) * (1 + 0.7 * Re^(1/5) * Pr^(1/3))
 *    + (1.33 - 2.4*ep + 1.2*ep^2) * Re^(7/10) * Pr^(1/3)
 *
 * For ep = 1 the correlation reduces to
 * Nu = 2 + 1.4 * Re^(1/5) * Pr^(1/3) + 0.13 * Re^(7/10) * Pr^(1/3).
 *
 * The Reynolds number is based on the superficial slip velocity, i.e.
 * Re = ep * rho_f * |u_f - u_p| * d_p / mu_f.
 *
 * The correlation is reported for 0.35 <= ep <= 1 and Re <= 10^5. No
 * clipping is applied inside this class: outside the validity range the
 * polynomials are simply evaluated as they are.
 *
 * Reference: D. J. Gunn, Transfer of heat or mass to particles in fixed
 * and fluidised beds, Int. J. Heat Mass Transfer, 21 (1978) 467-476.
 */
class Gunn
{
public:

    //- Type info

        TypeInfoNV("Gunn");

    //- constructors

        /**
         * @brief Constructor initializing the model from the simulation dict.
         * @param dict Sub-dictionary containing specific model coefficients.
         */
        explicit Gunn(const Foam::dictionary& dict);

        ~Gunn() = default;

    //- public methods

        /**
         * @brief Calculates the dimensionless Nusselt number.
         * @param Re Particle Reynolds number [-].
         * @param Pr Fluid Prandtl number [-].
         * @param ep Fluid volume fraction (porosity) [-].
         * @return The computed Nusselt number [-].
         */
        inline
        Foam::scalar dimlessNusselt(
            Foam::scalar Re,
            Foam::scalar Pr,
            Foam::scalar ep) const
        {
            const Foam::scalar Pr13 = Foam::pow(Pr, 1.0/3.0);

            return
                (7.0 - 10.0*ep + 5.0*ep*ep)
              * (1.0 + 0.7 * Foam::pow(Re, 0.2) * Pr13)
              + (1.33 - 2.4*ep + 1.2*ep*ep)
              * Foam::pow(Re, 0.7) * Pr13;
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

#endif // pFlow_Gunn_hpp
