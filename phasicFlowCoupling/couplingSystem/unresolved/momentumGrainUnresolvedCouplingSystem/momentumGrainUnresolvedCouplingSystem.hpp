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
 * @class momentumGrainUnresolvedCouplingSystem
 * @ingroup couplingSystem
 * 
 * @brief Momentum coupling system for spherical grains (course-grained particles) in unresolved particle-fluid flows.
 *
 * The `momentumGrainUnresolvedCouplingSystem` class implements Eulerian-Lagrangian momentum
 * coupling for **spherical grain particles** in unresolved simulations. This class extends
 * the unresolved coupling framework to handle arbitrarily shaped particles by incorporating
 * a **course-grain factor** that accounts for shape effects on drag and lift forces.
 *
 * ## Key Features
 *
 * - **Non-Spherical Particles**: Designed for grains and irregular particle shapes
 * - **Course-Grain Factor**: Applies particle-specific correction factor to force models
 * - **Shape-Dependent Drag**: Drag coefficients are adjusted based on particle course-grain factors
 * - **Momentum Coupling**: Calculates drag, lift, and virtual mass forces for non-spherical particles
 * - **Porosity Calculation**: Computes local fluid volume fraction in cells
 * - **Distribution Methods**: Leverages parent class methods to map grain properties to cells
 *
 * ## Physical Coupling Mechanisms
 *
 * ### Shape Correction via Course-Grain Factor
 * The course-grain factor $f_{cg}$ modifies the effective particle diameter for force calculations:
 * - Adjusts drag coefficient: $\beta_{adjusted} = \beta \cdot f_{cg}$
 * - Allows drag models (DiFelice, ErgunWenYu, Beetstra, etc.) to account for shape effects
 * - Specific to each particle (stored in `courseGrainFactor_` field)
 *
 * ### Momentum Source Term
 * Similar to spherical coupling but with shape-corrected forces:
 * - $S_p = \frac{1}{V_c} \sum_p f_{cg,p} \cdot w_{p,c} \frac{V_p \beta}{1-\alpha}$
 * - $S_u = -\frac{1}{V_c} \sum_p f_{cg,p} \cdot w_{p,c} \frac{V_p \beta}{1-\alpha} \mathbf{v}_p + \text{corrections}$
 *
 *

 *
 * ## Configuration Example
 *
 * ```C++
 * unresolved
 * {
 *     distribution  adaptiveGaussian;  
 *     adaptiveGaussianInfo { ... }
 *     
 *     porosity
 *     {
 *         method    distribution;
 *         alphaMin  0.2;
 *     }
 *     
 *     drag
 *     {
 *         model     DiFelice;          // Will be modulated by courseGrainFactor
 *         residualRe 1e-6;
 *     }
 *     
 *     lift
 *     {
 *         // not implemented yet for course-grain particles
 *     }
 * }
 * ```
 *
 * ## Supported Options
 *
 * - **Drag Models**: DiFelice, ErgunWenYu, Beetstra, Rong, Cello (with shape correction)
 * - **Lift Models**: not implemented for course-grain particles
 * - **Torque Models**: not implemented for course-grain particles
 * - **Heat Coupling**: not implemented (see notImplementedFunction calls)
 * - **Mass Coupling**: not implemented (see notImplementedFunction calls)
 * - **Distribution Methods**: All methods supported, adaptiveGaussian recommended for accuracy
 *
 *
 * @see unresolvedCouplingSystem, momentumSphereUnresolvedCouplingSystem
 * @see porosity, momentumInteraction
 *
 */

#ifndef __momentumGrainUnresolvedCouplingSystem_hpp__
#define __momentumGrainUnresolvedCouplingSystem_hpp__


#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "momentumInteraction.hpp"

namespace pFlow::coupling
{

class momentumGrainUnresolvedCouplingSystem
:
    public unresolvedCouplingSystem
{
private:

	/// Porosity calculator for computing local fluid volume fraction.
    uniquePtr<porosity>					porosity_ = nullptr;

	/// Drag, lift, and virtual mass force calculator for momentum coupling.
    momentumInteraction                 momentumInteraction_;

	/// Per-particle shape correction factors for spherical grain drag coefficient modification.
    Plus::realProcCMField		        courseGrainFactor_;

	/// Timer for monitoring porosity calculation performance.
    Timer 		porosityTimer_;

	/// Flag indicating if coupling requires distribution weights for field mapping.
    bool 		requiresDistribution_ = false;

protected:

	/// Override to distribute particle course-grain factors across MPI processors along with standard fields.
    bool distributeParticleFields() override;

public:

	/// Type information 
    TypeInfo("grainUnresolvedCouplingSystem<momentum>");

	/// Constructor initializing momentum coupling for spherical grain particles.
    momentumGrainUnresolvedCouplingSystem(
        word shapeTypeName,
        word couplingSystemType, 
        Foam::fvMesh& mesh,
        int argc, 
        char* argv[]);

	/// Deleted copy constructor to prevent unintended copying of coupling system.
    momentumGrainUnresolvedCouplingSystem(const momentumGrainUnresolvedCouplingSystem&) = delete;

	/// Deleted copy assignment operator to prevent unintended copying of coupling system.
    momentumGrainUnresolvedCouplingSystem& operator=(const momentumGrainUnresolvedCouplingSystem&) = delete;

	/// Deleted move constructor to prevent unintended moving of coupling system.
    momentumGrainUnresolvedCouplingSystem(momentumGrainUnresolvedCouplingSystem&&) = delete;

	/// Deleted move assignment operator to prevent unintended moving of coupling system.
    momentumGrainUnresolvedCouplingSystem& operator=(momentumGrainUnresolvedCouplingSystem&&) = delete;

	/// Virtual destructor for proper cleanup of derived classes.
    ~momentumGrainUnresolvedCouplingSystem() override = default;

	/// Virtual constructor macro for registering the derived class in the factory selector.
    add_vCtor
    (
        unresolvedCouplingSystem,
        momentumGrainUnresolvedCouplingSystem,
        word
    );

	/// Get const reference to per-particle course-grain factors.
	/// @return Const reference to realProcCMField containing course-grain factors for all particles.
    const Plus::realProcCMField& courseGrainFactor()const
    {
        return courseGrainFactor_;
    }

	/// Get mutable reference to per-particle course-grain factors.
	/// @return Mutable reference to realProcCMField for updating course-grain factors.
    Plus::realProcCMField& courseGrainFactor()
    {
        return courseGrainFactor_;
    }

	/// Calculate local porosity (fluid volume fraction) in each cell.
    void calculatePorosity() override;

	/// Calculate momentum coupling forces (drag, lift, virtual mass) between particles and fluid.
    void calculateMomentumCoupling() override;

	/// Calculate heat coupling (not yet implemented for grain particles).
    void calculateHeatCoupling() override;

	/// Calculate mass coupling (not yet implemented for grain particles).
    void calculateMassCoupling() override;

	/// Get the fluid volume fraction field (porosity) in cells.
	/// @return Const reference to the volScalarField representing fluid volume fraction.
    const Foam::volScalarField& alpha()const override
    {
        return porosity_().alpha();
    }

	/// Get the implicit coefficient of momentum source term (Sp in Sp*U+Su).
	/// @return Temporary volScalarField containing the implicit momentum source coefficient.
    Foam::tmp<Foam::volScalarField> Sp()const override;

	/// Get the explicit part of momentum source term (Su in Sp*U+Su).
	/// @return Temporary volVectorField containing the explicit momentum source.
    Foam::tmp<Foam::volVectorField> Su()const override;
    
	/// Get the heat source term (not yet implemented).
	/// @return Temporary volScalarField (returns null pointer; function not implemented).
    Foam::tmp<Foam::volScalarField> heatSource() const override
    {
        notImplementedFunction;
        return Foam::tmp<Foam::volScalarField>(nullptr);
    }

	/// Get the mass source term for species transport (not yet implemented).
	/// @param specieName Name of the chemical species or phase.
	/// @return Temporary volScalarField (returns null pointer; function not implemented).
    Foam::tmp<Foam::volScalarField> massSource(const word& specieName)const  override
    {
        notImplementedFunction;
        return Foam::tmp<Foam::volScalarField>(nullptr);
    }

	/// Get the particle shape type name.
	/// @return Word "grain" identifying the particle shape type (used for grain coupling).
    word shapeTypeName() const override
    {
        return "grain";
    }

	/// Get the coupling system type name.
	/// @return Word "momentum" identifying the coupling physics type.
    word couplingSystemType()const override
    {
        return "momentum";
    }

	/// Check if coupling requires cell-level distribution of coupling terms.
	/// @return True if distribution weights are needed; false if using cell centroid method.
    bool requireCellDistribution()const override
    {
        return requiresDistribution_;
    }
    
};


} // pFlow::coupling


#endif //__momentumGrainUnresolvedCouplingSystem_hpp__

