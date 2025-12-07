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
 * @class momentumSphereUnresolvedCouplingSystem
 * @ingroup couplingSystem
 * 
 * @brief Momentum coupling system for spherical particles in unresolved particle-fluid flows.
 *
 * The `momentumSphereUnresolvedCouplingSystem` class implements Eulerian-Lagrangian momentum
 * coupling for **spherical particles** in unresolved simulations. This class handles the exchange
 * of momentum between particles and fluid through drag, lift, and virtual mass forces, without
 * explicitly resolving the particle boundary layer in the fluid mesh.
 *
 * ## Key Features
 *
 * - **Spherical Particles**: Specialized for spherical particle geometries
 * - **Momentum Coupling**: Calculates and applies drag, lift, and virtual mass forces
 * - **Porosity Calculation**: Computes local fluid volume fraction (porosity) in cells
 * - **Multiple Drag Models**: Supports DiFelice, ErgunWenYu, Beetstra, Rong, and other models
 * - **Lift Force Models**: Saffman, Loth2008, Shi2019, and surface rotation torque models
 * - **Distribution Methods**: Leverages parent class distribution methods for smooth field mapping
 *
 * ## Physical Coupling Mechanisms
 *
 * ### Drag Force
 * Represents resistance to relative motion between particles and fluid:
 * - $F_D = \frac{V_p \beta}{1-\alpha}(\mathbf{U}-\mathbf{v}_p)$
 * - Uses dimensionless drag correlation: $\beta = \frac{18\mu_f \alpha(1-\alpha)}{d_p^2} \hat{f}^d(\alpha, Re)$
 *
 * ### Lift Force
 * Shear-induced and spin-induced lift effects:
 * - $\mathbf{F}_L = \frac{\pi d_p^2}{8} \rho_f |\mathbf{U}_{rel}|^2 (C_{L,shear} \mathbf{e}_{shear} + C_{L,spin} \mathbf{e}_{spin})$
 *
 * ### Virtual Mass Force
 * Represents inertia of displaced fluid (added mass effect)
 *
 * ## Momentum Source Term
 *
 * The fluid momentum equation receives source terms in the form:
 * - $S = S_p \mathbf{U} + S_u$
 * - Where $S_p$ is the implicit coefficient and $S_u$ is the explicit source
 * - These are computed from weighted contributions of all particles in each cell:
 *   - $S_p = \frac{1}{V_c} \sum_p w_{p,c} \frac{V_p \beta}{1-\alpha}$
 *   - $S_u = -\frac{1}{V_c} \sum_p w_{p,c} \frac{V_p \beta}{1-\alpha} \mathbf{v}_p + \text{lift, virtual mass}$
 *
 *
 * ## Configuration Example
 *
 * ```C++
 * unresolved
 * {
 *     distribution  Gaussian;      // or PCM, diffusion, adaptiveGaussian, etc.
 *     GaussianInfo  { ... }
 *     
 *     porosity
 *     {
 *         method    distribution;
 *         alphaMin  0.2;
 *     }
 *     
 *     drag
 *     {
 *         model     DiFelice;
 *         residualRe 1e-6;
 *     }
 *     
 *     lift
 *     {
 *         model     Loth2008;
 * 
 *         surfaceRotationTorque Shi2019;
 * 
 *         residualRe 1e-6;
 *     }
 * }
 * ```
 *
 * ## Supported Options
 *
 * - **Drag Models**: DiFelice, ErgunWenYu, Beetstra, Rong, Cello
 * - **Lift Models**: none, Saffman, Loth2008, Shi2019
 * - **surfaceRotationTorque models**: none, Shi2019, Loth2008, lowReynolds
 *
 * @see unresolvedCouplingSystem, momentumGrainUnresolvedCouplingSystem
 * @see porosity, momentumInteraction
 * @see drag, lift models in interaction subdirectory
 *
 * @note Heat and mass coupling are not yet implemented (see notImplementedFunction calls)
 */

#ifndef __momentumSphereUnresolvedCouplingSystem_hpp__
#define __momentumSphereUnresolvedCouplingSystem_hpp__


#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "momentumInteraction.hpp"

namespace pFlow::coupling
{

class momentumSphereUnresolvedCouplingSystem
:
	public unresolvedCouplingSystem
{
private:

	/// Porosity calculator for computing local fluid volume fraction.
	uniquePtr<porosity>					porosity_ = nullptr;

	/// Drag, lift, and virtual mass force calculator for momentum coupling.
    momentumInteraction                 momentumInteraction_;

	/// Timer for monitoring porosity calculation performance.
	Timer 		porosityTimer_;

	/// Flag indicating if coupling requires distribution weights for field mapping.
    bool 		requiresDistribution_ = false;

public:

	/// Type information 
    TypeInfo("sphereUnresolvedCouplingSystem<momentum>");

	/// Constructor initializing momentum coupling for spherical particles.
    momentumSphereUnresolvedCouplingSystem(
        word shapeTypeName,
        word couplingSystemType, 
        Foam::fvMesh& mesh,
        int argc, 
        char* argv[]);

	/// Deleted copy constructor to prevent unintended copying of coupling system.
    momentumSphereUnresolvedCouplingSystem(const momentumSphereUnresolvedCouplingSystem&) = delete;

	/// Deleted copy assignment operator to prevent unintended copying of coupling system.
    momentumSphereUnresolvedCouplingSystem& operator=(const momentumSphereUnresolvedCouplingSystem&) = delete;

	/// Deleted move constructor to prevent unintended moving of coupling system.
    momentumSphereUnresolvedCouplingSystem(momentumSphereUnresolvedCouplingSystem&&) = delete;

	/// Deleted move assignment operator to prevent unintended moving of coupling system.
    momentumSphereUnresolvedCouplingSystem& operator=(momentumSphereUnresolvedCouplingSystem&&) = delete;

	/// Virtual destructor for proper cleanup of derived classes.
    ~momentumSphereUnresolvedCouplingSystem() override = default;

	/// Virtual constructor macro for registering the derived class in the factory selector.
    add_vCtor
    (
        unresolvedCouplingSystem,
        momentumSphereUnresolvedCouplingSystem,
        word
    );

	/// Calculate local porosity (fluid volume fraction) in each cell.
    void calculatePorosity() override;

	/// Calculate momentum coupling forces (drag, lift, virtual mass) between particles and fluid.
    void calculateMomentumCoupling() override;

	/// Calculate heat coupling (not yet implemented for momentum coupling).
    void calculateHeatCoupling() override;

	/// Calculate mass coupling (not yet implemented for momentum coupling).
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
	/// @return Word "sphere" identifying the particle shape type.
    word shapeTypeName() const override
    {
        return "sphere";
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


#endif //__momentumSphereUnresolvedCouplingSystem_hpp__

