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

#ifndef pFlow_sphereReactionTransfer_hpp
#define pFlow_sphereReactionTransfer_hpp

#include <vector>

#include "reactionTransfer.hpp"
#include "distributionBase.hpp"

namespace pFlow::coupling
{

/**
 * @brief Lagrangian-to-Eulerian mapping of the reaction sources for spheres.
 *
 * Loops over the active particles and distributes their multi-species mass
 * sources and sinks into the OpenFOAM volumetric fields with OpenMP, using
 * the configured spatial distribution.
 */
class sphereReactionTransfer : public reactionTransfer
{
public:

    TypeInfo("sphereReactionTransfer<sphere>");

    //- constructors

        sphereReactionTransfer
        (
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty
        );

        ~sphereReactionTransfer() override = default;

        add_vCtor
        (
            reactionTransfer,
            sphereReactionTransfer,
            couplingSystem
        );

    //- public methods

        /**
         * @brief Distributes the reaction sources of the particles to the
         * CFD mesh cells.
         *
         * @param gasMassSource    [p * N_gas + g] = explicit mass source
         * [kg/s] of particle p, gas species g.
         * @param gasMassSourceSp  [p * N_gas + g] = implicit sink
         * coefficient [kg.m3/(kmol.s)].
         * @param cellDistribution Distribution operator (PCM or Gaussian).
         */
        void calculateReactionTransfer
        (
            const std::vector<real>&     gasMassSource,
            const std::vector<real>&     gasMassSourceSp,
            const distributionBase&      cellDistribution
        ) override;
};

} // namespace pFlow::coupling

#endif // pFlow_sphereReactionTransfer_hpp
