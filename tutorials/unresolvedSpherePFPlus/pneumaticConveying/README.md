# Pnumatic Conveying

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: unresolvedSpherePFPlus

## 0. Problem Definition

In this tutorial, we will simulate a gas-solid, horizontal pnumatic conveying system. Typically the length of the conduit in a real system is several meters. In a few meters from the injection point, the flow reaches a fully developed/well-stablished state. However, simulating such a long conduit with CFD-DEM coupling is computationally expensive. Therefore, in this tutorial, we will simulate a shorter conduit of length 1 m to demonstrate the stable operation of the pneumatic conveying in the fully developed region. To do so, we use a periodic/cyclic boundary condition in the flow direction (x-direction) to mimic a longer conduit. The siulation is perfored using the unresolved solver `unresolvedSpherePFPlus`.

<div align="center">
<b>
<img src="./pnueaticConveying.gif" alt="horizontal pnumatic conveying"/>
</b>
<b>

An animtation of the simulated horizontal pnumatic conveyying.
</b></div>

The length of the conduit is 1 m, with a square cross-section of 0.1x0.1 m². The gas (air) enters the conduit from the left side (inlet) and exits from the right side (outlet). The top and bottom walls of the conduit are no-slip walls. Solid particles with a density of 800 kg/m³ are injected into the conduit from a point source located 0.175 m downstream of the inlet with a rate of 40000 particles/s for 0.5 seconds (in total 20000 particles are inserted). Two sizes of particles are injected: 5 mm and 4 mm.

The simulation lasts for 1.5 seconds. In period between 0 s to 0.5 s, particles are inserted into the conduit and the rest of simulation continues without particle insertion. The gas flow under constant pressure drop (1500 Pa) across the tube is established from the start of the simulation. 

**A note on the pressure boundary condition**

In the CFD side of the simulaiton, the type of inlet and outlet boundaries are set to cyclic. Since, we cannot set an explicit value for the inlet velocity nor the pressure in the outlet, we cannot maintain the gas velocity in the simulation and the gas velocity usually decay to zero. To solve this problem, we have two options in OpenFOAM:

- Applying a pressure drop across the conduit using `uniformJump` boundary condition for pressure.

- Using `meanVelocityForce` option in fvOptions to apply a force in the flow direction to maintain the mean velocity of the gas at a specified value.

Here, we used the first option. In this way, we kept the presure difference across the conduit ends at a constat value and the average air velocity is ajusted accroding to the total flow resistance in the gas-solid system in the conduit. Obviousely, with increasing the particle concenteration in the system, the flow resistance increases and the average gas velocity decreases. At the operating conditions of this simulation (when the pressure drop across the conduit is around 1500 Pa and particle volume fraction of 10%), the average superficial gas velocity at stable condition is around 16.2 m/s.  

***

## 1. Performing the Simulation Using the Allrun Script

The `Allrun` script is designed to automate the simulation process using the `unresolvedSpherePFPlus` solver. It manages all essential steps, including mesh generation, CFD-DEM coupling simulation, and result conversion. 


## 2. Important parts

We must apply periodict/cyclic boundary conditions to both ends of the concute. In DEM part, this is possible through the definition of the domain. In `domainDict` file, the global domain is specified and it fitst the boundries conduit. since the periodic boundary is in x-direction, the type of boundaries in x-dicrection (`left` and `right`), is set to `periodic`.

<div align="center">
in <b>settings/domainDict</b>
</div>

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

In the CFD part, the type of inlet and outlet boundary patches should be set to `cyclic`, in the `blockMeshDict` file.

<div align="center">
in <b>system/blockMeshDict</b>
</div>

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

And in the `p` field, the following boundary conditions are applied to `inlet` and `outlet` patches: `uniformJump`. The value of jump is specified via `jumpTable`. Since the pressure in the solver is kinematic pressure and the value of density is 1.2 kg/m3, the value is set to 1250 m2/s2, which is equivalent to 1500 Pa.

<div align="center">
in <b>fluid/p</b>
</div>

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
