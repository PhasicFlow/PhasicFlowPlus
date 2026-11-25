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

#ifndef pFlow_massTransfer_hpp
#define pFlow_massTransfer_hpp

#include <vector>

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "porosity.hpp"
#include "procCMFields.hpp"

namespace pFlow
{

namespace coupling
{

class unresolvedCouplingSystem;
class fluidAveraging;
class solidAveraging;

/**
 * @brief Abstract base class for the gas-film mass transfer coefficient of
 * the particles.
 *
 * The coefficient k_c [m/s] of each particle is sent to the DEM, where it
 * enters the gas-film term of the particle model of the reactions.
 */
class massTransfer
{
private:

    //- private members

        /// Reference to the local void fraction (porosity) field.
        const porosity&     porosity_;

public:

    //- Type info

        TypeInfo("massTransfer");

        create_vCtor(
            massTransfer,
            couplingSystem,
            (
                const unresolvedCouplingSystem& uCS,
                const porosity&                 prsty
            ),
            (uCS, prsty)
        );

    //- constructors

        massTransfer(
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty);

        virtual ~massTransfer() = default;

    //- public methods

        inline
        const porosity& Porosity() const
        {
            return porosity_;
        }

        inline
        const auto& parCellIndex() const
        {
            return porosity_.parCellIndex();
        }

        inline
        const Foam::fvMesh& mesh() const
        {
            return porosity_.mesh();
        }

        inline
        const Foam::volScalarField& alpha() const
        {
            return static_cast<const Foam::volScalarField&>(porosity_);
        }

        const Foam::dictionary& dict() const;

        /// Calculates k_c [m/s] of each particle; 0 for a particle outside
        /// the fluid mesh of this processor.
        virtual void calculateMassTransferCoeff(
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    diameter,
            std::vector<real>&              massTransferCoeff) = 0;

    //- static members

        static
        uniquePtr<massTransfer> create(
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty);

        /// True when couplingProperties has the massTransferModel block.
        static
        bool found(const unresolvedCouplingSystem& uCS);

        static
        const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);

}; // massTransfer

} // coupling

} // pFlow

#endif // pFlow_massTransfer_hpp
