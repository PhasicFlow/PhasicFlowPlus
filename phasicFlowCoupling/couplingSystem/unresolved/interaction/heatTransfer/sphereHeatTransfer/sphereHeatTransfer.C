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

#include "sphereHeatTransfer.hpp"
#include <cmath> 

namespace pFlow
{
namespace coupling
{

//----------------------------- constructors ----------------------------------

template<typename NusseltClosureType>
sphereHeatTransfer<NusseltClosureType>::sphereHeatTransfer(
    const unresolvedCouplingSystem& uCS, 
    const porosity&                 prsty)
: 
    heatTransfer(uCS, prsty), 
    nusseltClosure_(this->dict()) 
{
    // Functor initialized with the sub-dictionary
}

//---------------------------- public methods ---------------------------------

template<typename NusseltClosureType>
void sphereHeatTransfer<NusseltClosureType>::calculateHeatTransfer(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    diameter,
    const Plus::realProcCMField&    parTemperature,
    const distributionBase&         cellDistribution,  
    Plus::realProcCMField&          Qp,                    
    const Plus::realProcCMField&    emissivity,            
    const Plus::realProcCMField&    radSumTemp,            
    const Plus::uint32ProcCMField&  radNumPrt,             
    Plus::realProcCMField&          QpRad)
{
    // Prevent source term accumulation from previous time steps
    this->setSuSpToZero();

    const auto& parCellInd = this->parCellIndex();
    const auto& Vcells     = this->mesh().V();
    
    // Retrieve required Eulerian fields from OpenFOAM's objectRegistry
    const auto& nu = 
        this->mesh().template lookupObject<Foam::volScalarField>("nu");
        
    const auto& rho = 
        this->mesh().template lookupObject<Foam::volScalarField>("rho");
        
    const auto& Cp = 
        this->mesh().template lookupObject<Foam::volScalarField>("Cp");
        
    const auto& kappa = 
        this->mesh().template lookupObject<Foam::volScalarField>("kappa");
        
    const auto& Tfld = 
        this->mesh().template lookupObject<Foam::volScalarField>("T");
    
    size_t numPar     = parCellInd.size();
    const auto& alpha = this->alpha();
    
    auto fluidVel = fluidVelocity.fieldSpan();
    auto solidVel = parVelocity.fieldSpan();

    auto& Su = this->Su();
    auto& Sp = this->Sp();

    // ------------------------------------------------------------------------
    // Compile-time constants
    // ------------------------------------------------------------------------
    constexpr Foam::scalar stefanBoltzmann = 5.670374419e-8;
    constexpr Foam::scalar pi = Foam::constant::mathematical::pi;

    // Memory & Size Validation
    const bool hasRadiation = 
        (radNumPrt.size() == numPar) && 
        (radSumTemp.size() == numPar) && 
        (emissivity.size() == numPar);

    // ------------------------------------------------------------------------
    // Unresolved Assumption Guard (Executed only once per simulation)
    //
    // Checks the LARGEST particle diameter across the ENTIRE domain (all
    // processors), not just the first particle in this processor's local
    // array. The "numPar > 0" guard is intentionally NOT used here: both
    // the local-maximum loop and the collective reduction below must run
    // identically on every processor on every call, regardless of how many
    // particles happen to be mapped to this processor right now.
    // ------------------------------------------------------------------------
    if (!this->unresolvedWarningIssued_)
    {
        Foam::scalar globalMaxDiameter = 0.0;

        for (size_t parIndx = 0; parIndx < numPar; parIndx++)
        {
            globalMaxDiameter = 
                Foam::max(globalMaxDiameter, diameter[parIndx]);
        }

        // Collective, all-processor reduction: combine every processor's
        // local maximum into a single domain-wide maximum diameter.
        Foam::reduce(globalMaxDiameter, Foam::maxOp<Foam::scalar>());

        const Foam::scalar minV        = Foam::gMin(Vcells);
        const Foam::scalar minCellSize = Foam::pow(minV, 1.0/3.0);

        if (globalMaxDiameter > minCellSize)
        {
            WarningInFunction
                << "Largest particle diameter in the domain (" 
                << globalMaxDiameter
                << ") exceeds the minimum fluid cell size (" << minCellSize
                << "). The unresolved CFD-DEM assumption (V_cell >> V_p) "
                << "may be violated!" << Foam::endl;
        }

        this->unresolvedWarningIssued_ = true;
    }

    // ------------------------------------------------------------------------
    // Thermodynamic Bounds
    // Limits fetched outside the loop to avoid thread contention.
    //
    // alphaMin is read directly from the porosity model instead of
    // a second, independent dictionary entry, so this bound can never
    // silently disagree with the one used to compute the alpha field.
    // ------------------------------------------------------------------------
    const Foam::scalar kappaMin = 
        this->dict().lookupOrDefault<Foam::scalar>("kappaMin", 1e-15);
        
    const Foam::scalar muMin = 
        this->dict().lookupOrDefault<Foam::scalar>("muMin", 1e-15);
        
    const Foam::scalar dpMin = 
        this->dict().lookupOrDefault<Foam::scalar>("dpMin", 1e-15);
        
    const Foam::scalar alphaMin = this->Porosity().alphaMin();
    
    const Foam::scalar Tmin = 
        this->dict().lookupOrDefault<Foam::scalar>("Tmin", 1.0);

    // ------------------------------------------------------------------------
    // High-Performance OpenMP Loop
    // Thread-safety relies on atomic operations inside cellDistribution
    // ------------------------------------------------------------------------
    #pragma omp parallel for schedule(dynamic)
    for (size_t parIndx = 0; parIndx < numPar; parIndx++)
    {
        auto cellIndx = parCellInd[parIndx];
        
        if (cellIndx < 0) continue;

        // Extract and clamp required scalars locally
        const Foam::scalar k_val = Foam::max(kappa[cellIndx], kappaMin);
        
        const Foam::scalar mui = 
            Foam::max(nu[cellIndx] * rho[cellIndx], muMin);
            
        const Foam::scalar dp = Foam::max(diameter[parIndx], dpMin);
        
        const Foam::scalar ef = 
            Foam::min(Foam::max(alpha[cellIndx], alphaMin), 1.0); 
            
        const Foam::scalar Tf = Foam::max(Tfld[cellIndx], Tmin);
        
        const Foam::scalar Tp_val = Foam::max(parTemperature[parIndx], Tmin);

        // Calculate relative slip velocity
        const Foam::vector up{
            solidVel[parIndx].x(), 
            solidVel[parIndx].y(), 
            solidVel[parIndx].z()
        };
        const Foam::vector ur = fluidVel[parIndx] - up;
        
        // Dimensionless Numbers
        const Foam::scalar Re = 
            ef * rho[cellIndx] * Foam::mag(ur) * dp / mui;
            
        const Foam::scalar Pr = (Cp[cellIndx] * mui) / k_val;
        
        // Convection Calculation
        Foam::scalar Nu = nusseltClosure_(Re, Pr, ef);

        // NOTE: this floor is intentionally left disabled. Enforcing
        // Nu >= 2 unconditionally would break the "none" heat-transfer
        // closure (noneHeatTransfer), whose entire purpose is to return
        // Nu = 0 so that convective coupling can be switched off for
        // isolated testing of other mechanisms (conduction, PFP,
        // radiation). Real closures (e.g. RanzMarshall) already satisfy
        // Nu >= 2 on their own and do not need this floor.
        // Nu = Foam::max(Nu, 2.0);

        const Foam::scalar hAp_conv = Nu * k_val * pi * dp;
        
        // Heat transfer rate injected INTO the particle (Explicit)
        Qp[parIndx] = hAp_conv * (Tf - Tp_val);

        // Particle-to-Fluid source terms (Implicit formulation for stability)
        const Foam::scalar spConv = -hAp_conv;                
        const Foam::scalar suConv =  hAp_conv * Tp_val;       

        // --------------------------------------------------------------------
        // Radiation Calculation
        //
        // Modelling assumption: the carrier fluid is radiatively transparent.
        // It neither emits nor absorbs radiation. Radiative energy travels
        // directly between solid particles; the neighbourhood-averaged
        // temperature T_p_avg (blended with the local fluid temperature Tf
        // only as a proxy for the unresolved "background" beyond the nearest
        // few neighbours) drives ONLY the particle's own absorbed radiative
        // flux QpRad. This must NEVER be added to spConv/suConv, since those
        // feed the fluid energy equation (Su/Sp) and the fluid must have
        // zero net energy exchange through this mechanism.
        //
        // A particle with no radiating neighbours (radNumPrt == 0) correctly
        // exchanges no radiative heat at all: there is nothing to radiate to
        // or receive from, and the transparent fluid cannot substitute for a
        // missing neighbour.
        // --------------------------------------------------------------------
        if (hasRadiation && radNumPrt[parIndx] > 0)
        {
            const Foam::scalar nPrt = 
                static_cast<Foam::scalar>(radNumPrt[parIndx]);
                
            const Foam::scalar T_p_avg = 
                Foam::max(radSumTemp[parIndx] / nPrt, 1.0);

            const Foam::scalar T_env = ef * Tf + (1.0 - ef) * T_p_avg;

            const Foam::scalar emis = 
                Foam::min(Foam::max(emissivity[parIndx], 0.0), 1.0);
                
            const Foam::scalar Ap = pi * dp * dp;

            const Foam::scalar hAp_rad = emis * stefanBoltzmann * Ap 
                                       * (T_env*T_env + Tp_val*Tp_val) 
                                       * (T_env + Tp_val);

            if (QpRad.size() > parIndx)
            {
                QpRad[parIndx] = hAp_rad * (T_env - Tp_val);
            }

            // Deliberately NOT added to spConv/suConv
        }

        // Distribute ONLY the convective exchange to the fluid mesh.
        // Radiative exchange (QpRad, above) never reaches the fluid
        // energy equation - it is carried to the DEM side instead.
        cellDistribution.distributeValue_OMP(
            parIndx, 
            cellIndx, 
            Su, 
            suConv);
            
        cellDistribution.distributeValue_OMP(
            parIndx, 
            cellIndx, 
            Sp, 
            spConv);
    }

    // Convert total heat exchange to proper volumetric source terms [W/m^3]
    forAll(Vcells, i)
    {
        Su[i] /= Vcells[i];
        Sp[i] /= Vcells[i];
    }

    // Apply Gaussian or Cell-based smoothing
    cellDistribution.smoothenField(Sp);
    cellDistribution.smoothenField(Su);

    // Synchronize across processor boundaries
    Sp.correctBoundaryConditions();
    Su.correctBoundaryConditions();
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling
} // pFlow



