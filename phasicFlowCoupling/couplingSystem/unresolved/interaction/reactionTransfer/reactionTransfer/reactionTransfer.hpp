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

#ifndef pFlow_reactionTransfer_hpp
#define pFlow_reactionTransfer_hpp

#include <vector>

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "processorPlus.hpp"
#include "porosity.hpp"
#include "procCMFields.hpp"
#include "distributionBase.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

/**
 * @brief Abstract base class for mapping the DEM reaction sources onto the
 * CFD mesh. Holds one volumetric mass source (Su) and one implicit sink
 * coefficient (Sp) field per gas species.
 */
class reactionTransfer
{
private:

    //- private members

        /// Explicit mass source of each gas species [kg/(m3.s)].
        Foam::PtrList<Foam::volScalarField> speciesSu_;

        /// Implicit sink coefficient of each gas species [kg/(kmol.s)].
        Foam::PtrList<Foam::volScalarField> speciesSp_;

        /// Reference to the local void fraction (porosity) field.
        const porosity& porosity_;

        /// Gas species names (productGas of couplingProperties).
        Foam::wordList gasSpeciesNames_;

protected:

    //- protected methods

        /// Resets the explicit and implicit source fields to zero.
        void setSuSpToZero();

        /// Non-const access to the explicit source fields.
        inline Foam::PtrList<Foam::volScalarField>& speciesSu()
        {
            return speciesSu_;
        }

        /// Non-const access to the implicit source fields.
        inline Foam::PtrList<Foam::volScalarField>& speciesSp()
        {
            return speciesSp_;
        }

public:

    TypeInfo("reactionTransfer");

    //- constructors

        reactionTransfer
        (
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty
        );

        virtual ~reactionTransfer() = default;

        create_vCtor
        (
            reactionTransfer,
            couplingSystem,
            (
                const unresolvedCouplingSystem& uCS,
                const porosity&                 prsty
            ),
            (uCS, prsty)
        );

    //- public methods

        inline const porosity& Porosity() const { return porosity_; }

        inline const auto& parCellIndex() const
        {
            return porosity_.parCellIndex();
        }

        inline const Foam::fvMesh& mesh() const { return porosity_.mesh(); }

        inline const Foam::wordList& gasSpeciesNames() const
        {
            return gasSpeciesNames_;
        }

        inline const Foam::PtrList<Foam::volScalarField>& speciesSu() const
        {
            return speciesSu_;
        }

        inline const Foam::PtrList<Foam::volScalarField>& speciesSp() const
        {
            return speciesSp_;
        }

        const Foam::dictionary& dict() const;

        /// Distributes the particle reaction sources onto the mesh fields.
        virtual void calculateReactionTransfer
        (
            const std::vector<real>&     gasMassSource,
            const std::vector<real>&     gasMassSourceSp,
            const distributionBase&      cellDistribution
        ) = 0;

        static uniquePtr<reactionTransfer> create
        (
            const unresolvedCouplingSystem& uCS,
            const porosity&                 prsty
        );

        static const Foam::dictionary& getDict
        (
            const unresolvedCouplingSystem& uCS
        );
};

} // namespace pFlow::coupling

#endif // pFlow_reactionTransfer_hpp
