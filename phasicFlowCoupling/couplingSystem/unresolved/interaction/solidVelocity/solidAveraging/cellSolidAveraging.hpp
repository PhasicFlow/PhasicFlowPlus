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

#ifndef __cellSolidAveraging_hpp__
#define __cellSolidAveraging_hpp__

#include "solidAveraging.hpp"

namespace pFlow::coupling
{

template<typename DistributionType>
class cellSolidAveraging
:
    public solidAveraging<DistributionType>
{
public:

    using SolidAveragingType = solidAveraging<DistributionType>;

    using CellSolidAveragingType = cellSolidAveraging<DistributionType>

private:

    Foam::volVectorField        cellAvField_;

    Plus::realx3ProcCMField     calcParAvField_;

public:

    // type info
    TypeInfoTemplate111("solidAveraging", DistributionType, "cell");

    cellSolidAveraging(
        const word& type,
        const couplingMesh& cMesh,
        Plus::realx3ProcCMField& particleField);

    ~cellSolidAveraging()override =default;

    add_vCtor
	(
		SolidAveragingType,
		CellSolidAveragingType,
		word	
	);
    
    /// @brief perform mass-based averaging (if applicable)
    /// @param distributor the cellDistribution system
    /// @param prsty porosity object  
    void calculate(
        const Plus::realx3ProcCMField& particleField
        const DistributorType& cellDistributor,
        const porosity& prsty) override;
    
    /// @brief perform number-based averaging (if applicable)
    /// @param distributor the cellDistribution system
    /// @param prsty porosity object 
    void calculateNumberBased(
        const Plus::realx3ProcCMField& particleField
        const DistributorType& cellDistributor,
        const porosity& prsty) override;
    
    bool requireCellDistribution()const override
    {
        return false;
    }

};

}

#include "cellSolidAveraging.cpp"

#endif //__cellSolidAveraging_hpp__