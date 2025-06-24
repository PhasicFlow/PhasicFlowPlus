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

#ifndef __particleMapping_hpp__
#define __particleMapping_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

#include "scatteredCommunicationPlus.hpp"
#include "centerMassField.hpp"
#include "box.hpp"

namespace pFlow::Plus
{
class procDEMSystem;
}

namespace pFlow::coupling
{

class couplingMesh;


class particleMapping
:
    public 	Plus::procCommunication
{
    /// center point of particles in this processor 
    Plus::centerMassField 		centerMass_;

    /// The mesh domain expansion ratio
    /// Number of times the largest particle in the simulation
    Foam::scalar		domainExpansionRatio_;

    /// Time intervals between domain updates
    Foam::scalar 		domainUpdateInterval_;

    /// Last time of domain update
    Foam::scalar 		lastTimeUpdated_ = 0;

    Plus::scatteredCommunication<real> 		realScatteredComm_;

	Plus::scatteredCommunication<realx3> 	realx3ScatteredComm_;

    Plus::scatteredCommunication<uint32>    uint32ScatteredComm_;

    /// box containing the mesh for all processors
    Plus::procVector<box> meshBoxes_;

    /// If everything is constructed for the first time
    bool                firstConstructed_ = false;


public:

    particleMapping(
        const Foam::dictionary& dict);

    /// Domain extension ratio
    inline
    auto domainExpansionRatio()const
    {
        return domainExpansionRatio_;
    }

    /// Check if the domain should be updated at time t
    /// In fluid loop, the current time is dt ahead of coupling time, 
    /// so the function is notified to consider this. 
    bool checkForDomainUpdate
    (
        Foam::scalar t, 
        Foam::scalar fluidDt
    )const;

    bool update(
        Foam::scalar t,
        Foam::scalar fluidDt,
        Plus::procDEMSystem& pDEMSystem,
        const couplingMesh& cMesh);

    inline
    Plus::centerMassField& centerMass()
    {
        return centerMass_;
    }

    inline
    const Plus::centerMassField& centerMass()const
    {
        return centerMass_;
    }

    Plus::scatteredCommunication<real>& realScatteredComm() 
    {
        return realScatteredComm_;
    }

    Plus::scatteredCommunication<realx3>& realx3ScatteredComm() 
    {
        return realx3ScatteredComm_;
    }

    Plus::scatteredCommunication<uint32>& uint32ScatteredComm() 
    {
        return uint32ScatteredComm_;
    }

    inline 
	const Plus::procVector<box>& meshBoxes()const
	{
		return meshBoxes_;
	}
};

}

#endif //__particleMapping_hpp__
