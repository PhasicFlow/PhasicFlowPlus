/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
-----------------------------------------------------------------------------*/

#ifndef pFlow_heatTransfer_hpp
#define pFlow_heatTransfer_hpp

#include "virtualConstructor.hpp"
#include "processorPlus.hpp"
#include "porosity.hpp"
#include "procCMFields.hpp"

namespace pFlow
{
namespace coupling
{

class distributionBase;
class unresolvedCouplingSystem;
class fluidAveraging;
class solidAveraging;

/**
 * @class heatTransfer
 * @brief Abstract base class for all interphase heat transfer models.
 *
 * @details
 * This class defines the strict interface that any specific heat transfer 
 * closure model (e.g., Ranz-Marshall, Gunn) must implement. It manages 
 * the allocation and access to the volumetric source terms (Su and Sp) that 
 * will be injected directly into the continuous phase (OpenFOAM) energy matrix.
 */
class heatTransfer
{
public:

    //- Type info

        TypeInfo("heatTransfer");

        create_vCtor(
            heatTransfer,
            couplingSystem,
            (
                const unresolvedCouplingSystem& uCS, 
                const porosity&                 prsty
            ),
            (uCS, prsty)
        );

private:

    //- private members

        /// @brief Explicit heat source field (RHS of the energy equation) 
        /// [W/m^3].
        Foam::volScalarField    heatSu_;
        
        /// @brief Implicit heat source coefficient (Diagonal of energy matrix) 
        /// [W/(m^3.K)].
        Foam::volScalarField    heatSp_;
        
        /// @brief Reference to the local void fraction (porosity) field.
        const porosity&         porosity_;

protected:

    //- protected methods

        /// @brief Resets both explicit (Su) and implicit (Sp) fields to zero.
        void setSuSpToZero();

        /// @brief Non-const access to explicit source field.
        inline
        Foam::volScalarField& Su()
        {
            return heatSu_;
        }
        
        /// @brief Non-const access to implicit source field.
        inline
        Foam::volScalarField& Sp()
        {
            return heatSp_;
        }

public:

    //- constructors

        heatTransfer(
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty);

        virtual ~heatTransfer() = default;

    //- public methods

        inline
        const porosity& Porosity() const
        {
            return porosity_;
        }
        
        inline
        const auto& parCellIndex() const
        {
            return porosity_.parCellIndex();
        }
        
        inline
        const Foam::fvMesh& mesh() const
        {
            return porosity_.mesh();
        }
        
        inline
        const auto& cMesh() const
        {
            return porosity_.cMesh();
        }
        
        inline
        const Foam::volScalarField& alpha() const 
        { 
            return static_cast<const Foam::volScalarField&>(porosity_); 
        }
        
        inline
        const Foam::volScalarField& Su() const
        {
            return heatSu_;
        }
        
        inline
        const Foam::volScalarField& Sp() const
        {
            return heatSp_;
        }

        const Foam::dictionary& dict() const;

        /**
         * @brief Core calculation method for specific closure models.
         * @param fluidVelocity Interpolated fluid velocity at particles.
         * @param parVelocity   Particle velocity array.
         * @param diameter      Particle diameter array.
         * @param parTemperature Particle temperature array.
         * @param cellDistribution Mapping object for distribution.
         * @param Qp            [OUT] Computed convective heat rate array.
         * @param emissivity    Particle emissivity array.
         * @param radSumTemp    Sum of neighboring temperatures for radiation.
         * @param radNumPrt     Number of neighbors considered in radiation.
         * @param QpRad         [OUT] Computed radiative heat rate array.
         */
        virtual void calculateHeatTransfer(
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    diameter,
            const Plus::realProcCMField&    parTemperature,
            const distributionBase&         cellDistribution,  
            Plus::realProcCMField&          Qp,
            const Plus::realProcCMField&    emissivity,
            const Plus::realProcCMField&    radSumTemp,
            const Plus::uint32ProcCMField&  radNumPrt,
            Plus::realProcCMField&          QpRad) = 0;

        static const Foam::dictionary& getDict(
            const unresolvedCouplingSystem& uCS);    
        
        static uniquePtr<heatTransfer> create(
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty);

}; 

} // coupling
} // pFlow

#endif // pFlow_heatTransfer_hpp
