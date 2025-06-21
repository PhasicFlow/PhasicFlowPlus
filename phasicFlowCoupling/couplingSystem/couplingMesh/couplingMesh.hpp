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
#include "OFCompatibleHeader.hpp"
#include "indexedOctree.H"
#include "treeDataCell.H"



// phasicFlow
#include "procCMField.hpp"
#include "box.hpp"


namespace pFlow::coupling
{

class couplingMesh
{
private:

	// - private data members

		/// Reference to fvMesh
		Foam::fvMesh& 					mesh_;

		/// number of cells in the mesh
		Foam::label 					nCells_;

		/// cell indices of particles in this processor 
		Plus::procCMField<Foam::label> 	parCellIndex_;

		/// number of particles found in this mesh
		int32 							numInMesh_ = 0;

		/// Octree for cell search
		mutable uniquePtr<Foam::indexedOctree<Foam::treeDataCell>>
							cellTreeSearch_ = nullptr;

		/// Actual bounding box around mesh points
		mutable box 		meshBox_;	

		/// cell decomposition mode
		Foam::polyMesh::cellDecomposition 		cellDecompositionMode_;

	// - member functions

		/// Calculate the actual bounding box mesh based on points
		void calculateBox()const;

		/// Reset the search tree structure
		void resetTree()const;

		void mapParticles();

public:

	// - Constructors

		couplingMesh(
			const Foam::dictionary& dict,
			Foam::fvMesh& mesh,
			const Plus::centerMassField& centerMass
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
        void update();

        /// cell index of each particle center 
        inline
		const Plus::procCMField<Foam::label>& parCellIndex()const
        {
            return parCellIndex_;
        }

		inline
		const Plus::centerMassField& centerMass()const
		{
			return parCellIndex_.centerMass();
		}

        int32 numInMesh()const
        {
            return numInMesh_;
        }

        /// Report (output) number of center mass points found in all processors 
		/// It is effective only in master processor 
		void reportNumInMesh()const;
		
        /// check if a point is located in a specific cell
		bool pointInCell(
			const Foam::point& p, 
			Foam::label celli)const;

		bool pointSphereInCell(
			const Foam::point& p, 
			Foam::scalar rad, 
			Foam::label celli, 
			bool& sphereInCell) const;

		bool pointSphereInCell(
			const Foam::point& p, 
			Foam::scalar smallRad, 
			Foam::scalar largeRad, 
			Foam::label celli,
			bool& smallInCell, 
			bool& largeInCell)const;


		Foam::label findPointInCellTree(
			const Foam::point& p, 
			Foam::label cellId)const;

		Foam::label findPointSphereInCellTree(
			const Foam::point& p,
			Foam::scalar rad,
			Foam::label cellId,
			bool& sphereInCell)const;

		Foam::label findPointSphereInCellTree(
			const Foam::point& p,
			Foam::scalar radSmall,
			Foam::scalar radLarge,
			Foam::label cellId,
			bool& smallInCell,
			bool& largeInCell)const;

		Foam::label
		findCellTree(const realx3& p, Foam::label cellId)const;

		template<unsigned Size>
		void findPointsInCells(
			const Foam::FixedList<Foam::point, Size>& points, 
			Foam::label cntrCellId, 
			Foam::label& nCellIds,
			Foam::FixedList<Foam::label, Size>& cellIds)const
		{
			nCellIds = 0;
			for(auto i=0u; i<Size; i++)
			{
				if(auto id = findPointInCellTree(points[i], cntrCellId); cntrCellId!=id && id !=-1)
				{
					cellIds[nCellIds] = id;
					nCellIds++;
				}
			}
		}
		

		template<unsigned Size>
		void findPointsInCells(
			const Foam::FixedList<realx3, Size>& points, 
			Foam::label cntrCellId, 
			Foam::label& nCellIds,
			Foam::FixedList<Foam::label, Size>& cellIds)const
		{
			nCellIds = 0;
			for(auto i=0u; i<Size; i++)
			{
				if(auto id = findCellTree(points[i], cntrCellId); cntrCellId!=id && id !=-1)
				{
					cellIds[nCellIds] = id;
					nCellIds++;
				}
			}
		}


		Foam::labelList findSphere(Foam::label cellId, Foam::scalar radius)const;
		

};

}



#endif //__couplingMesh_hpp__



/*template<unsigned Size>
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
		}*/