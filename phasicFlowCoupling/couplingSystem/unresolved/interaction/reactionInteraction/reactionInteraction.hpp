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

#ifndef pFlow_reactionInteraction_hpp
#define pFlow_reactionInteraction_hpp

#include <utility>
#include <vector>

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "Timer.hpp"
#include "procCMFields.hpp"
#include "reactionTransfer.hpp"
#include "PCM.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;
class porosity;

/**
 * @brief Applies the discrete reaction sources of the DEM particles to the
 * CFD mesh.
 */
class reactionInteraction
{
private:

    /// Void fraction (porosity) coupling.
    const porosity&             porosity_;

    /// True for "reactionExchange distribution", false for "cell".
    bool                        reactionExchangeDistribute_;

    /// True when the distribution weights of the particles are needed.
    bool                        requireCellDistribution_ = false;

    /// Mapping of the particle sources onto the mesh fields.
    uniquePtr<reactionTransfer> reactionTransfer_;

    /// Single-cell (PCM) mapping, used when no distribution is requested.
    uniquePtr<PCM>              noDistribution_ = nullptr;

    /// Timer of the reaction interaction.
    Timer                       reactionInteractionTimer_;

public:

    TypeInfo("reactionInteraction");

    reactionInteraction
    (
        const unresolvedCouplingSystem& uCS,
        const porosity& prsty
    );

    virtual ~reactionInteraction() = default;

    /// Explicit gas mass source of each species [kg/(m3.s)], for YEqn.
    inline const Foam::PtrList<Foam::volScalarField>& speciesSu() const
    {
        return std::as_const(*reactionTransfer_).speciesSu();
    }

    /// Implicit sink coefficient of each species [kg/(kmol.s)], for YEqn.
    inline const Foam::PtrList<Foam::volScalarField>& speciesSp() const
    {
        return std::as_const(*reactionTransfer_).speciesSp();
    }

    /// Names of the gas species, in the order of speciesSu() and speciesSp().
    inline const Foam::wordList& gasSpeciesNames() const
    {
        return reactionTransfer_->gasSpeciesNames();
    }

    inline const porosity& Porosity() const { return porosity_; }

    const unresolvedCouplingSystem& uCS() const;

    inline bool requireCellDistribution() const
    {
        return requireCellDistribution_;
    }

    const Foam::dictionary& dict() const;

    static const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);

    virtual void calculateCoupling
    (
        const std::vector<real>&     gasMassSource,
        const std::vector<real>&     gasMassSourceSp
    );
};

} // pFlow::coupling

#endif // pFlow_reactionInteraction_hpp
