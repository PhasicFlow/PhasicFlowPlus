# Pnumatic Conveying

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: unresolvedSpherePFPlus

## 0. Problem definition

This tutorial simulates a gas–solid horizontal pneumatic conveying system.

In real pneumatic conveying lines, the conduit length is typically several meters. After a short entrance region downstream of the injection point, the flow reaches a fully developed (quasi-steady) state. Simulating a full-length line with CFD–DEM coupling is computationally expensive, so here we simulate a shorter conduit (1 m) that represents a fully developed section.

To mimic an effectively long conduit, we apply periodic/cyclic boundary conditions in the streamwise direction (the $x$ direction). The simulation is performed with the unresolved solver `unresolvedSpherePFPlus`.

<div align="center">
<b>
<img src="./pnueaticConveying.gif" alt="horizontal pnumatic conveying"/>
</b>
<b>

An animtation of the simulated horizontal pnumatic conveyying.
</b></div>

The conduit length is 1 m, with a square cross-section of $0.1\times0.1\ \mathrm{m}^2$. The gas (air) enters from the left side (inlet) and exits from the right side (outlet). The top and bottom walls are no-slip walls.

Solid particles (density $800\ \mathrm{kg/m^3}$) are injected from a point source located 0.175 m downstream of the inlet. The insertion rate is 40,000 particles/s for 0.5 s (total of 20,000 particles). Two particle sizes are injected: 4 mm and 5 mm diameter spheres.

The simulation lasts for 2 s. Particles are inserted from 0 s to 0.5 s, and the simulation continues afterwards without further insertion. A constant pressure drop (1500 Pa) is applied across the periodic tube from the start of the simulation.

### A note on the pressure boundary condition

On the CFD side, the inlet and outlet patches are set to `cyclic`. With purely cyclic boundaries, you cannot prescribe an inlet velocity or an outlet pressure directly, so the mean flow can decay toward zero unless a driving force is applied. In OpenFOAM, there are two common ways to drive the flow:

- Apply a pressure drop across the conduit using the `uniformJump` boundary condition for pressure.

- Use `meanVelocityForce` in `fvOptions` to apply a body force in the flow direction and maintain a specified mean gas velocity.

Here we use the first option. With `uniformJump`, the pressure difference between the two periodic ends is held constant, and the mean gas velocity adjusts according to the total flow resistance of the gas–solid mixture in the conduit.

As the particle concentration increases, the flow resistance increases and the mean gas velocity decreases. For the operating conditions of this case (pressure drop $\approx 1500\ \mathrm{Pa}$ and particle volume fraction $\approx 10\%$), the superficial gas velocity at steady operation is about 15.8 m/s.

***

## 1. Running the case

Run the case using the provided script `runThisCase`. It automates the full workflow:

- mesh generation (`blockMesh`)
- DEM setup (`particlesPhasicFlow`, `geometryPhasicFlow`)
- copying initial fluid fields into `0/`
- coupled solve (`unresolvedSpherePFPlus`)
- post-processing export (`pFlowToVTK`)


## 2. Key settings

We apply periodic/cyclic boundary conditions at both ends of the conduit.

On the DEM side, this is done by defining a domain that matches the conduit and setting the streamwise boundaries to periodic. In `settings/domainDict`, the global domain matches the conduit bounds, and since periodicity is in the $x$ direction, the `left` and `right` boundaries are set to `periodic`.

In **`settings/domainDict`**:

```C++
// Simulation domain: every particle that goes outside this domain will be deleted
globalBox
{
    min ( 0.0 0.0 0.0);
    max ( 1.0 0.1 0.1);
}

boundaries
{
    neighborListUpdateInterval 10;
    updateInterval 1;

    left // negative x
    {
        type periodic; 
    }
    right // positive x
    {
        type periodic; 
    }
    bottom // negative y
    {
        type exit; 
    }
    top // positive y
    {
        type exit; 
    }
    rear // negative z
    {
        type exit; 
    }
    front // positive z
    {
        type exit; 
    }
}
```

On the CFD side, the corresponding patches must be `cyclic` in `system/blockMeshDict`.

In **`system/blockMeshDict`**:

```C++
// some lines are missing 

boundary
(
    // some lines are missing

    top
    {
        type wall;
        //some lines are missin
    }
    bottom 
    {
        type wall;
        // some lines are missing
    }
    outlet
    {
        type cyclic;
        neighbourPatch inlet;
        faces
        (
            (2 6 5 1)
        );
    }
    inlet
    {
        type cyclic;
        neighbourPatch outlet;
        faces
        (
            (0 4 7 3)
        );
    }
    front
    {
        type wall;
        // some lines are missing
    }
    rear
    {
        type wall;
        // some lines are missing 
    }
);
```

Finally, in the `p` field, the `inlet` and `outlet` patches use `uniformJump`. The jump magnitude is set via `jumpTable`.

Because the solver uses kinematic pressure and the fluid density is $\rho=1.2\ \mathrm{kg/m^3}$, a pressure drop of 1500 Pa corresponds to a kinematic pressure jump of $\Delta p/\rho = 1500/1.2 \approx 1250\ \mathrm{m^2/s^2}$.

In **`fluid/p`**:

```C++
dimensions      [ 0 2 -2 0 0 0 0 ];

internalField   uniform 0; 

boundaryField
{
    bottom
    {
        type            fixedFluxPressure;
    }
    top
    {
        type            fixedFluxPressure;
    }
    rear
    {
        type            fixedFluxPressure;
    }
    front
    {
        type            fixedFluxPressure;
    }
    outlet
    {
        type            uniformJump;
        patchType       cyclic;
        jumpTable       constant 1250; // ~ to 1500 Pa (density of fluid is 1.2 kg/m3) 
        relax           0.5;
        value           uniform 0;    
    }
    inlet
    {
        type            uniformJump;
        patchType       cyclic;
    }
}
```
