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

#ifndef pFlow_sphereMassTransfer_hpp
#define pFlow_sphereMassTransfer_hpp

#include "massTransfer.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow
{

namespace coupling
{

/**
 * @brief Gas-film mass transfer coefficient of spherical particles.
 *
 * The Sherwood number comes from one of the Nusselt closures of
 * heatTransferModel, used with the Schmidt number in place of the Prandtl
 * number (heat and mass transfer analogy):
 *
 *     Sh  = closure(Re, Sc, ep)
 *     k_c = Sh * D / d_p,   D = nu / Sc
 *
 * Re = ep * |u_f - u_p| * d_p / nu is the Reynolds number of
 * sphereHeatTransfer; Sc is read from transportProperties, where YEqn takes
 * the same value for the gas diffusivity.
 */
template<typename SherwoodClosureType>
class sphereMassTransfer
:
    public massTransfer
{
public:

    /// Alias of the template type for the registration macro.
    using SphereMassTransferType = sphereMassTransfer<SherwoodClosureType>;

private:

    //- private members

        /// Closure (functor) for the Sherwood number.
        SherwoodClosureType     sherwoodClosure_;

public:

    //- Type info

        TypeInfoTemplate11("sphereMassTransfer", SherwoodClosureType);

    //- constructors

        sphereMassTransfer(
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty);

        ~sphereMassTransfer() override = default;

        add_vCtor(
            massTransfer,
            SphereMassTransferType,
            couplingSystem);

    //- public methods

        void calculateMassTransferCoeff(
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    diameter,
            std::vector<real>&              massTransferCoeff) override;

}; // sphereMassTransfer

} // coupling

} // pFlow

// Template implementation
#include "sphereMassTransfer.C"

#endif // pFlow_sphereMassTransfer_hpp
