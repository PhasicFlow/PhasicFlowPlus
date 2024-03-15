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

#ifndef __couplingMesh_hpp__ 
#define __couplingMesh_hpp__

// from OpneFOAM
#include "fvMesh.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

#include <map>

// phasicFlow
#include "box.hpp"


namespace pFlow::coupling
{

class couplingMesh
{
protected:

	// - protected data members

		/// Reference to fvMesh
		Foam::fvMesh& 		mesh_;

		/// Is mesh static
		bool 				meshStatic_;

		/// Cctree for cell search
		mutable uniquePtr<Foam::indexedOctree<Foam::treeDataCell>>
							cellTreeSearch_ = nullptr;

		/// The mesh domain expansion ratio
		/// Number of times the largest particle in the simulation
		Foam::scalar		domainExpansionRatio_;

		/// Time intervals between domain updates
		Foam::scalar 		domainUpdateInterval_;

		/// Last time of domain update
		Foam::scalar 		lastTimeUpdated_ = 0;

		/// Mesh decomposition mode for locating points in the mesh
		/// Options are: 
		/// 	1) facePlanes
		/// 	2) cellTets
		/// 	3) faceDiagonalTriangles
		/// 	4) faceCenterTriangles
		Foam::word 			decompositionMode_;
		
		/// Actual bounding box around mesh points
		mutable box 		meshBox_;

		/// cell decomposition mode
		Foam::polyMesh::cellDecomposition 		cellDecompositionMode_;

		/// If everything is constructed for the first time
		bool 				firstConstruction_ = false;
		

	// - protected member functions

		/// Calculate the actual bounding box mesh based on points
		void calculateBox()const;

		/// Reset the search tree structure
		void resetTree()const;


	/*Foam::label findCellSeed(
		const Foam::point& loc,
    	const Foam::label seedCellId);*/

public:

	// - Constructors

		couplingMesh(
			const Foam::dictionary& dict,
			Foam::fvMesh& mesh
			);

		couplingMesh(const couplingMesh&)=delete;

		couplingMesh(couplingMesh&&)=delete;

		couplingMesh& operator = (const couplingMesh&)=delete;

		couplingMesh& operator = (couplingMesh&&)=delete;		

		~couplingMesh()=default;


	// - Methods 

		/// Retrun const ref to mesh
		inline
		const auto& mesh()const
		{
			return mesh_;
		}

		/// Domain extension ratio
		inline
		auto domainExpansionRatio()const
		{
			return domainExpansionRatio_;
		}

		/// Return cosnt ref to mesh box
		inline
		const auto& meshBox()const
		{
			return meshBox_;
		}
		
		/// Is mesh dynamic?
		inline 
		auto dynamic()const
		{
			return mesh_.dynamic();
		}
		
		/// Is mesh moving?
		inline 
		auto moving()const
		{
			return mesh_.moving();
		}
		
		/// Is topology of the mesh is changing?
		inline 
		auto topoChanging()const
		{
			return mesh_.topoChanging();
		}
		
		/// Is mesh changing
		inline
		auto changing()const
		{
			return mesh_.changing();
		}
		
		/// Update the coupling mesh components
		/// This should always be called before any mesh querries
		/// and after mesh motion (if any).
		void update(Foam::scalar t, Foam::scalar fluidDt);

		/// Check if the domain should be updated at time t
		/// In fluid loop, the current time is dt ahead of coupling time, 
		/// so the function is notified to consider this. 
		bool checkForDomainUpdate
		(
			Foam::scalar t, 
			Foam::scalar fluidDt, 
			bool insideFluidLoop = true
		);

		Foam::label
		findCellTree(const realx3& p, Foam::label cellId)const;

		template<unsigned Size>
		void findPointsInCells(
			const Foam::FixedList<realx3, Size>& points, 
			Foam::label cntrCellId, 
			Foam::label& nCellIds,
			Foam::FixedList<Foam::label, Size>& cellIds)
		{
			nCellIds = 0;
			for(auto i=0; i<Size; i++)
			{
				if(auto id = findCellTree(points[i], cntrCellId); cntrCellId!=id && id !=-1)
				{
					cellIds[nCellIds] = id;
					nCellIds++;
				}
			}
		}

		template<unsigned Size>
		void findPointsInCellsMap(
			const Foam::FixedList<realx3, Size>& points, 
			Foam::label cntrCellId, 
			Foam::label& nCellIds,
			std::map<Foam::label, Foam::label>& cellIds)
		{
			nCellIds = 0;
			cellIds.clear();
			for(auto i=0; i<Size; i++)
			{
				if(auto id = findCellTree(points[i], cntrCellId); cntrCellId!=id && id !=-1)
				{
					if( auto cIdIter = cellIds.find(id); cIdIter != cellIds.end())
					{
						cIdIter->second++;
					}
					else
					{
						cellIds.insert({id,1});
					}
					
					
					nCellIds++;
				}
			}
		}

};

}



#endif //__couplingMesh_hpp__
