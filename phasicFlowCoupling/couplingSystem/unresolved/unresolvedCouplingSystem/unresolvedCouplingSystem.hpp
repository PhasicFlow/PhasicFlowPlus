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
 * @class unresolvedCouplingSystem
 * @ingroup couplingSystem
 * 
 * @brief Base class for Eulerian-Lagrangian coupling systems in unresolved particle-fluid flows.
 *
 * The `unresolvedCouplingSystem` class serves as the abstract base for implementing unresolved
 * coupling between particles and fluid in multiphase flows. 
 *
 * ## Key Responsibilities
 *
 * - **Particle Distribution**: Manages how particle-scale properties (velocity, forces, etc.)
 *   are mapped to fluid cells using configurable distribution methods (e.g., Gaussian,
 *   diffusion, PCM, sub-division methods).
 * - **Mesh Management**: Handles the coupling mesh that bridges particles and fluid cells.
 * - **Abstract Interface**: Defines the interface for derived classes to implement specific
 *   coupling physics (momentum, heat, mass transfer).
 *
 * ## Architecture
 *
 * This class is part of a three-level hierarchy:
 * - **unresolvedCouplingSystem** (base): Defines distribution and mesh management
 * - **Derived Classes** (e.g., momentumSphereUnresolvedCouplingSystem, momentumGrainUnresolvedCouplingSystem):
 *   Implement specific physical coupling mechanisms
 *
 * ## Distribution Methods
 *
 * The distribution method determines how particle data is mapped to cells:
 * - **PCM** (Particle Centroid Method): Direct assignment to containing cell
 * - **Gaussian**: Gaussian kernel-based distribution
 * - **GaussianIntegral**: Volume-weighted Gaussian integration
 * - **adaptiveGaussian**: Adaptive kernel based on cell-to-particle size ratio
 * - **diffusion**: Laplacian diffusion-based smoothing
 * - **subDivision9/29**: Geometric sub-division methods
 *
 *
 * @see couplingSystem, momentumSphereUnresolvedCouplingSystem, momentumGrainUnresolvedCouplingSystem
 * @see distributionBase, porosity, momentumInteraction
 *
 */

#ifndef __unresolvedCouplingSystem_hpp__
#define __unresolvedCouplingSystem_hpp__


#include "couplingSystem.hpp"
#include "distributionBase.hpp"
#include "virtualConstructor.hpp"


namespace pFlow::coupling
{


class unresolvedCouplingSystem
:
    public couplingSystem
{
private:

    /// Distribution method for mapping particle properties to fluid cells.
    uniquePtr<distributionBase>         distribution_;

public:

    /// Constructor initializing the unresolved coupling system with particle and fluid mesh information.
    unresolvedCouplingSystem(
        word        shapeTypeName,
        word        couplingSystemType, 
        Foam::fvMesh& mesh,
        int         argc, 
        char*       argv[]);

    /// Deleted copy constructor to prevent unintended copying of coupling system.
    unresolvedCouplingSystem(const unresolvedCouplingSystem&) = delete;

    /// Deleted copy assignment operator to prevent unintended copying of coupling system.
    unresolvedCouplingSystem& operator=(const unresolvedCouplingSystem&) = delete;

    /// Deleted move constructor to prevent unintended moving of coupling system.
    unresolvedCouplingSystem(unresolvedCouplingSystem&&) = delete;

    /// Deleted move assignment operator to prevent unintended moving of coupling system.
    unresolvedCouplingSystem& operator=(unresolvedCouplingSystem&&) = delete;

    /// Virtual destructor for proper cleanup of derived classes.
    ~unresolvedCouplingSystem() override = default;

    /// Virtual constructor macro for registering derived classes in the factory selector.
    create_vCtor
    (
        unresolvedCouplingSystem,
        word,
        (
            word shapeTypeName,
            word couplingSystemType, 
            Foam::fvMesh& mesh,
            int argc, 
            char* argv[]
        ),
        (shapeTypeName, couplingSystemType, mesh, argc, argv)
    );

    /// Access the `unresolved` coupling configuration subdictionary.
    const Foam::dictionary& unresolvedDict()const;

    /// Const reference to the distribution method for mapping particle properties to cells.
    const distributionBase& distribution()const
    {
        return *distribution_;
    }

    /// Get the name of the distribution method (e.g., Gaussian, diffusion, PCM).
    const Foam::word distributionMethodName()const
    {
        return distribution_->distributionMethodName();
    }

    /// Const reference to the cell indices where particles are located onto cells.
    inline
    const Plus::procCMField<Foam::label>& parCellIndex()const
    {
        return cMesh().parCellIndex();
    }

    /// Update distribution weights based on current particle positions.
    void updateDistributionWeights()
    {
        if(distribution_)
            distribution_->updateWeights(this->particleDiameter());
    }

    /// Pure virtual method to calculate local fluid volume fraction (porosity) in cells.
    virtual
    void calculatePorosity() =0;

    /// Pure virtual method to calculate momentum coupling between particles and fluid.
    virtual
    void calculateMomentumCoupling() = 0;

    /// Pure virtual method to calculate heat coupling between particles and fluid.
    virtual 
    void calculateHeatCoupling() = 0;

    /// Pure virtual method to calculate mass coupling between particles and fluid.
    virtual
    void calculateMassCoupling() = 0;

    /// Pure virtual method returning the fluid volume fraction field (porosity).
    /// @return Const reference to the volScalarField alpha representing porosity.
    virtual 
    const Foam::volScalarField& alpha()const =0;

    /// Pure virtual method returning the implicit coefficient of momentum source term (Sp in Sp*U+Su).
    /// @return Temporary volScalarField containing the implicit momentum source coefficient.
    virtual
    Foam::tmp<Foam::volScalarField> Sp()const =0;

    /// Pure virtual method returning the explicit part of momentum source term (Su in Sp*U+Su).
    /// @return Temporary volVectorField containing the explicit momentum source.
    virtual
    Foam::tmp<Foam::volVectorField> Su()const =0;

    /// Pure virtual method returning the heat source term for the energy equation.
    /// @return Temporary volScalarField containing the heat source.
    virtual 
    Foam::tmp<Foam::volScalarField> heatSource() const = 0;

    /// Pure virtual method returning the mass source term for the species transport equation.
    /// @param specieName Name of the chemical species or phase.
    /// @return Temporary volScalarField containing the mass source for the specified specie.
    virtual
    Foam::tmp<Foam::volScalarField> massSource(const word& specieName)const = 0;

    /// Pure virtual method returning the particle shape type name (e.g., sphere, grain).
    /// @return Word identifier for the particle shape type.
    virtual 
    word shapeTypeName() const= 0;

    /// Pure virtual method returning the coupling system type (e.g., momentum, heatMomentum, massHeatMomentum).
    /// @return Word identifier for the coupling system type.
    virtual 
    word couplingSystemType()const = 0;

    /// Pure virtual method indicating if the coupling requires cell-level distribution of coupling terms.
    /// @return True if distribution weights are needed; false if using cell centroid method.
    virtual 
    bool requireCellDistribution()const = 0;

    /// Static factory method to create derived coupling system instances based on shape and coupling types.
    /// @param shapeTypeName Type of particle shape (e.g., sphere, grain).
    /// @param couplingSystemType Type of coupling (e.g., momentum).
    /// @param mesh OpenFOAM fluid mesh.
    /// @param argc Command line argument count.
    /// @param argv Command line arguments.
    /// @return Unique pointer to the created unresolved coupling system.
    static
    uniquePtr<unresolvedCouplingSystem> create
    (
        word shapeTypeName,
        word couplingSystemType, 
        Foam::fvMesh& mesh,
        int argc, 
        char* argv[]
    );

}; 

} // pFlow::coupling

#endif //__unresolvedCouplingSystem_hpp__
