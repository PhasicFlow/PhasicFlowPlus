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
#ifndef __porosity_hpp__ 
#define __porosity_hpp__

// from OpenFOAM
#include "dictionary.H"
#include "volFields.H"

// from phasicFlow
#include "virtualConstructor.hpp"


// from phasicFlow-coupling
#include "couplingMesh.hpp"
#include "procCMFields.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

/**
 * Interface class for porosity calculation models 
 * 
 * Interface class for calculating porosity of fluid in each fluid cell
 * based on the various models
 */
class porosity
:
	public Foam::volScalarField
{
protected:

	/// Minimum fluid porosity allowed in each cell
	Foam::scalar					alphaMin_;

	/// Reference to couplingMesh
	const couplingMesh& 			cMesh_;

	/// Reference to diameter of particles in this processor 
	const Plus::realProcCMField&  	particleDiameter_;

	/// cell indices of particles in this processor 
	Plus::procCMField<Foam::label> 	parCellIndex_;

	int32 							numInMesh_ = 0;

	void mapCenters();

	bool centersMappedBefore_ = false;


protected:

	void setAlphaMin(Foam::scalar newAlphaMin)
	{
		alphaMin_ = Foam::max(Foam::min(newAlphaMin, 1.0), 0.001);
	}	

	template<typename DistributorType>
	Foam::tmp<Foam::volScalarField::Internal>
		calculateSolidVol(const DistributorType& distributor)
	{
		auto solidVolTmp = Foam::volScalarField::Internal::New
		(
			"solidVol",
			mesh(),
		 	Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		);

		
		auto& solidVol = solidVolTmp.ref();
		const Foam::label numPar = centerMass().size();

		#pragma omp parallel for
		for(Foam::label i=0; i<numPar; i++)
		{
			Foam::scalar pVol = pFlow::Pi/6 *
					Foam::pow(particleDiameter_[i], static_cast<real>(3.0));

			const Foam::label cellId = parCellIndex_[i];
			if( cellId >= 0 )
			{
				distributor.distributeValue_OMP(i, cellId, solidVol, pVol);				
			}
		}

		return solidVolTmp;
	}



public:

	/// Type info
	TypeInfo("porosity");

	// - Construcotrs 

		/// Construct from dictionary
		porosity(
			const unresolvedCouplingSystem& CS,
			const couplingMesh& 			cMesh,
			const Plus::realProcCMField& 	parDiam);

		/// No copy
		porosity(const porosity&) = delete;

		/// No copy assignment
		porosity& operator=(const porosity&) = delete;

		/// No move
		porosity(porosity&&) = delete;

		/// No move assignment
		porosity& operator=(porosity&&) = delete;

		/// Destructor 
		virtual ~porosity() = default;

		/// Virtual constructor 
		create_vCtor
		(
			porosity,
			dictionary,
			(
				const unresolvedCouplingSystem& CS,
				const couplingMesh& 			cMesh,
				const Plus::realProcCMField& 	parDiam
			),
			(CS, cMesh, parDiam)
		);

	// - Methods
		
		/// Const access to alpha
		inline
		const Foam::volScalarField& alpha()const
		{
			return *this;
		}

		inline 
		const Foam::fvMesh& mesh()const
		{
			return cMesh_.mesh();
		}
		/// Return alphaMin
		inline
		real alphaMin()const
		{
			return alphaMin_;
		}

		/// Return cell index of particles 
		inline 
		const auto& particleCellIndex()const
		{
			return parCellIndex_;
		}

		/// Retrun center mass position of particles 
		inline
		const auto& centerMass()const
		{
			return particleDiameter_.centerMass();
		}

		/// Retrun diameter of particles
		inline 
		const auto& particleDiameter()const
		{
			return particleDiameter_;
		}

		void mapCentersBeforeCalcPorosity()
		{
			centersMappedBefore_ = true;
			mapCenters();
		}
		/// Calculate porosity based on particles positions 
		void calculatePorosity();

		/// Fill the internal field of alpha
		virtual
		bool internalFieldUpdate() = 0;

		/// Return number of center mass points found in this mesh (processor)
		int32 numInMesh()const
		{
			return numInMesh_;
		}

		/// Report (output) number of center mass points found in all processors 
		/// It is effective only in master processor 
		void reportNumInMesh()const;

		virtual 
		bool requireCellDistribution()const
		{
			return false;
		}

	/// Construct the class based on input in dict 
	static
	uniquePtr<porosity> create(
		const unresolvedCouplingSystem& CS,
		const couplingMesh& 			cMesh,
		const Plus::realProcCMField& 	parDiam);
	

}; 

} // pFlow::coupling


#endif
