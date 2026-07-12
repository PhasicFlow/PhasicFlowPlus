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

/**
 * @class drag
 * @brief Base class for drag force models in unresolved Eulerian-Lagrangian
 *        coupling.
 *
 * This abstract base class defines the interface for drag force models
 * used to compute momentum exchange between particles and fluid in
 * unresolved particle-fluid simulations.
 *
 * @details
 * This class manages the source terms (Su_, Sp_), explicit and implicit parts,
 * for the fluid momentum equation and provides the interface for
 * computing drag forces using various closure models (e.g., DiFelice,
 * ErgunWenYu, Beetstra).
 */

#ifndef __drag_hpp__ 
#define __drag_hpp__

// from OpenFOAM
#include "OFCompatibleHeader.hpp"

// from phasicFlow
#include "virtualConstructor.hpp"

// from phasicFlow-coupling
#include "processorPlus.hpp"
#include "porosity.hpp"



namespace pFlow::coupling
{

class distributionBase;
class unresolvedCouplingSystem;
class fluidAveraging;
class solidAveraging;

class drag
{
private:

    /// @brief Explicit source term vector for momentum equation.
    Foam::volVectorField       Su_;

    /// @brief Implicit source term coefficient for momentum equation.
    Foam::volScalarField       Sp_;

    /// @brief Reference to porosity object (void fraction).
	const porosity& 			porosity_;

    /// @brief Reference to pressure field.
	const Foam::volScalarField& p_;

    /// @brief Flag indicating compressible flow solver.
	bool 						isCompressible_ = false;

protected:

    /// @brief Initialize source terms to zero.
    void setSuSpToZero();

    /// @brief Access explicit source term (non-const).
    inline
    Foam::volVectorField& Su()
    {
        return Su_;
    }

    /// @brief Access implicit source term (non-const).
    inline
    Foam::volScalarField& Sp()
    {
        return Sp_;
    }

public:

    /// Type info
    TypeInfo("drag");

    /// @brief Constructor from unresolved coupling system and porosity.
    drag(
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty);

    /// @brief Destructor.
    virtual ~drag() = default;

    create_vCtor
    (
        drag,
        couplingSystem,
        (
            const unresolvedCouplingSystem& uCS, 
            const porosity& 				prsty
        ),
        (uCS, prsty)
    );

    /// @brief Compute pressure gradient contribution to drag.
    Foam::tmp<Foam::volVectorField> 
    pressureGradient(const Foam::volScalarField& rho)const;

    /// @brief Check if flow is in compressible regime.
    inline
    bool isCompressible()const
    {
        return isCompressible_;
    }

    /// @brief Get constant reference to porosity object.
    const porosity& Porosity()const
    {
        return porosity_;
    }

    /// @brief Get particle-to-cell index mapping.
    inline
    const auto& parCellIndex()const
    {
        return porosity_.parCellIndex();
    }

    /// @brief Get reference to finite volume mesh.
    inline
    const Foam::fvMesh& mesh()const
    {
        return porosity_.mesh();
    }

    /// @brief Get reference to coupling mesh.
    inline 
    const auto& cMesh()const
    {
        return porosity_.cMesh();
    }

    /// @brief Get porosity field as volume scalar field.
    const Foam::volScalarField& alpha()const
    {
        return static_cast<const Foam::volScalarField&>(porosity_);
    }

    /// @brief Get constant explicit source term vector.
    inline
    const Foam::volVectorField& Su()const
    {
      return Su_;
    }

    /// @brief Get constant implicit source term scalar.
    inline
    const Foam::volScalarField& Sp()const
    {
      return Sp_;
    }

    /// @brief Get reference to drag dictionary.
    const Foam::dictionary& dict()const;

    /// @brief Compute drag force and populate particle force field.
    virtual 
    void calculateDragForce(
        const fluidAveraging& 	        fluidVelocity,
        const solidAveraging& 	        parVelocity,
        const Plus::realProcCMField& 	diameter,
        const distributionBase&         cellDistribution,    
        Plus::realx3ProcCMField& 		particleForce) = 0;
    
    /// @brief Get drag configuration dictionary from coupling system.
    static
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);    
    
    /// @brief Factory method to create drag model instance.
    static
    uniquePtr<drag> create
    (
        const unresolvedCouplingSystem& uCS, 
        const porosity& 				prsty
    );

}; 

} // pFlow::coupling


#endif
