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

#ifndef __drag_hpp__ 
#define __drag_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "processorPlus.hpp"
#include "UnresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "self.hpp"



namespace pFlow::coupling
{

template<typename DistributorType>
class drag
{
private:

	const porosity& 		porosity_;

	const DistributorType&   cellDistribution_;

	const Foam::volScalarField& p_;

	Foam::volVectorField 		Su_;

	Foam::volScalarField 		Sp_;

	bool 						isCompressible_ = false;

protected:

	/// private member functions 
	void setSuSpToZero();

	inline
	Foam::volVectorField& Su()
	{
		return Su_;
	}

	inline
	Foam::volScalarField& Sp()
	{
		return Sp_;
	}

public:

	// type info
	TypeInfoTemplate11("drag", DistributorType);

	drag(
		const UnresolvedCouplingSystem<DistributorType>& uCS, 
		const porosity& 				prsty);

	virtual ~drag() = default;

	create_vCtor
	(
		drag,
		couplingSystem,
		(
			const UnresolvedCouplingSystem<DistributorType>& uCS, 
			const porosity& 				prsty
		),
		(uCS, prsty)
	);

	Foam::tmp<Foam::volVectorField> 
	pressureGradient(const Foam::volScalarField& rho)const;

	inline
	const Foam::volVectorField& Su()const
	{
		return Su_;
	}

	inline
	const Foam::volScalarField& Sp()const
	{
		return Sp_;
	}
	
	inline
	bool isCompressible()const
	{
		return isCompressible_;
	}
	
	const porosity& Porosity()const
	{
		return porosity_;
	}

	inline
	const auto& parCellIndex()const
	{
		return porosity_.parCellIndex();
	}

	inline
	const Foam::fvMesh& mesh()const
	{
		return porosity_.mesh();
	}

	inline 
	const auto& cMesh()const
	{
		return porosity_.cMesh();
	}

	const Foam::volScalarField& alpha()const
	{
		return static_cast<const Foam::volScalarField&>(porosity_);
	}

	const DistributorType& cellDistribution()const
	{
		return cellDistribution_;
	}

	virtual 
	void calculateDragForce(
		const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce) = 0;
	
	
	virtual 
	Foam::scalar 
		dimlessDrag(Foam::scalar Re, Foam::scalar ep) =0;

	void calculateDragForce1(const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce);

	virtual
	bool requireCellDistribution()const =0;
	

	static
	uniquePtr<drag> create
	(
		const UnresolvedCouplingSystem<DistributorType>& uCS, 
		const porosity& 				prsty
	);
		
	
}; 

} // pFlow::coupling

#include "drag.C"

#endif
