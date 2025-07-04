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
#ifndef __GaussianOld_hpp__ 
#define __GaussianOld_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "typeInfo.hpp"

// from phasicFlowPlus
#include "procCMField.hpp"


namespace pFlow::coupling
{

class couplingMesh;

class GaussianOld
{
private:
	
	// radius of circule for cell neighbor list 
	Foam::scalar 				neighborLength_;

	Foam::scalar 				lengthExtent_;

	Plus::procCMField< std::vector<Foam::scalar> > 		weights_;

	// list of nieghbor cells for each cell 
	Foam::labelListList  		neighborList_;
	
	// boundary cells 
	Foam::List<std::pair<Foam::label, Foam::label>> 	boundaryCell_;

	const Foam::fvMesh&			mesh_;

protected:
	
	void constructLists(const Foam::scalar searchLen);

	void parseNeighbors(
		const Foam::label 		targetCelli,
		const Foam::vector& 	targetCellCentre, 
		const Foam::scalar 		searchLen,
		const Foam::scalar 		celli,
		std::set<Foam::label>& 	finalList,
		const Foam::label 		layerNumber);

public:

	/// Type info
	TypeInfoNV("GaussianOld");

	/// Construc from dictionary 
	GaussianOld(
		Foam::dictionary 		dict, 
		const couplingMesh& 	cMesh,
		const Plus::centerMassField& centerMass);

	/// Destructor
	~GaussianOld() = default;

	inline
	auto neighborLength()const
	{
		return neighborLength_;
	}

	void updateWeights(
		const Plus::procCMField<Foam::label> & parCellIndex,
		const Plus::procCMField<real> & parDiameter);

	inline
	void distributeValue_OMP(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volScalarField::Internal& internalField,
		const Foam::scalar& val)const
	{
		const auto& weightsI = weights_[parIndx];
		const Foam::labelList& neighbors = neighborList_[parCellIndx];
		forAll(neighbors, j)
		{
			#pragma omp atomic
			internalField[neighbors[j]] += val*weightsI[j];
		}
	}

	inline
	void distributeValue_OMP(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volVectorField::Internal& internalField,
		const Foam::vector& val )const
	{
		const auto& weightsI = weights_[parIndx];
		const Foam::labelList& neighbors = neighborList_[parCellIndx];
		forAll(neighbors, j)
		{
			const auto v = weightsI[j]* val;
			auto& tv = internalField[neighbors[j]]; 
			#pragma omp atomic
			tv.x() += v.x();
			
			#pragma omp atomic
			tv.y() += v.y();

			#pragma omp atomic
			tv.z() += v.z();
		}
	}

	inline 
	void distributeValue(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volScalarField::Internal& internalField,
		const Foam::scalar& val)const
	{
		const auto& weightsI = weights_[parIndx];
		const Foam::labelList& neighbors = neighborList_[parCellIndx];
		forAll(neighbors, j)
		{
			internalField[neighbors[j]] += val*weightsI[j];
		}
	}

	inline
	void distributeValue(
		Foam::label parIndx, 
		Foam::label parCellIndx,
		Foam::volVectorField::Internal& internalField,
		const Foam::vector& val)const
	{
		const auto& weightsI = weights_[parIndx];
		const Foam::labelList& neighbors = neighborList_[parCellIndx];
		forAll(neighbors, j)
		{
			internalField[neighbors[j]] += weightsI[j]* val;
		}
	}

	void smoothenField(Foam::volScalarField& field)const
	{}

	void smoothenField(Foam::volVectorField& field)const
	{}

	bool requireCellDistribution()const 
	{
		return true;
	}
}; 

} // pFlow::coupling


#endif
