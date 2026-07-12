# Cylinderical Fluidized Bed

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: unresolvedSpherePFPlus

## 0. Problem Definition

In this tutorial, we will simulate a gas-solid fluidized bed using the unresolved solver `unresolvedSpherePFPlus`. The cylinderical fluidized bed dimensions are 0.1 m in diameter and 0.2 m in height. We are going to 50K 2-mm particles with density of 1000 kg/m³. Gas is uniformly injected from the bottom of the bed at a superficial velocity of 1.2 m/s. The simulation runs for a total of 10 seconds, with the 0.25 s dedicated to the initial packing of particles (pure DEM simulation) and the remaining 9.75 seconds to the fluidized bed simulation.

<div align="center">
<b>
<img src="./fluidbed-snapshot.jpeg" alt="Fluidized bed" style="width: 400px;"/>
</b>
<b>

A visualization of a gas-solid fluidized bed with the gas field colored based on velocity.
</b></div>

***

## 1. Performing the Simulation Using the Allrun Script

The `Allrun` script is designed to automate the simulation process for the gas-solid fluidized bed using the `unresolvedSpherePFPlus` solver. It manages all essential steps, including mesh generation, DEM simulation, CFD-DEM coupling simulation, and result conversion. The first phase of the simulation is dedicated to particle settling, which is a pure DEM simulation using the `sphereGranFlow` solver. Following this phase, the CFD-DEM simulation is performed for the remaining 9.75 seconds. To execute the simulation, follow these steps:

### Step 1: Execute the `Allrun` Script

1. Navigate to the `fluidizedbed` directory.
2. Run the `Allrun` script by executing the following command:

   ```sh
   ./Allrun
   ```

   This script automates the entire simulation workflow, including mesh generation, DEM simulation, CFD-DEM coupling, and result conversion.

### Step 2: Understand the Folder Structure

The simulation folder structure is divided into two main categories: folders containing files related to DEM parameters and folders containing files related to CFD and coupling parameters:

- **DEM-related folders**:
  - **`settings/`**: Contains configuration files for the DEM simulation.
  - **`caseSetup/`**: Includes files for setting up the simulation and physical parameters for particles.

- **CFD-DEM-related folders**:
  - **`FluidField/`**: Holds the initial field data (`alpha`, `p`, `U`) used for the CFD-DEM simulation.
  - **`constant/`**: Contains constant properties for the fluid and parameters for coupling (CFD-DEM).
  - **`system/`**: Contains files for setting up CFD simulation parameters.

### Step 3: Visualize the Results

Once the simulation is complete, the results are converted to VTK format for visualization. The VTK files can be found in the `./VTK` folder. To visualize the results, use the following command:

```sh
paraview foam.foam &
```

Open the `foam.foam` file in ParaView to view the CFD results. For DEM results, open the `./VTK/particles.vtk.series` file.

## 2. Description of setup files

To learn about how to set up a DEM simulation, please refer to the [tutorial page](https://github.com/PhasicFlow/phasicFlow/wiki/Tutorials) of PhasicFlow and other online documents along side this package. Also, you can refer to OpenFOAM tutorials to learn about how to set up a CFD simulation. 

Mesh generation for the cylinder is performed using `blockMesh` utility of OpenFOAM. Although it is possible to use third-party mesh generation tools and then convert the mesh into OpenFOAM mesh. But here, we still use `blockMesh` utility. In `blockMeshDict` file, you will find these few lines, through which you can control the parameters for mesh generation:

```C++
// * * * User INPUTS  * * * * * * * * * *

  // Gemoetry 
  radius      0.05;
  height      0.2;
  cellSize    0.003;
  innerArc    yes; // yes/no

  // Patch information 
  inletPatchName      inlet;
  outletPatchName     outlet;
  wallPatchName       cylinderWall;
// * * * * * * * * * * * * * * * * * * *

// the rest of blockMeshDict file is missing here ... 
```
The parameters `radius`, `height`, and `cellSize` are used to define radius, hight and approximate edge cell size in the final mesh. You can also change the name of the patchs in the final mesh trough keys `inputPatchName`, `outletPatchName`, and `wallPatchName`. 

For geometry generation in DEM side (after mesh generation), we can use the utility that converts OpenFOAM patch surfaces into actual DEM walls. In the `settings/geometryDict` file, just set the type of the surface to `foamPatchWall`. This utility reads OpenFOAM mesh and converts the patches to DEM surfaces.  

```C++
surfaces
{
    inlet
    {
        type        foamPatchWall;   // type of the wall

        patch       inlet;      
        
        material    wallMaterial;   // material name of this wall
    }

    outlet
    {
        type        foamPatchWall;   // type of the wall

        patch       outlet;         
        
        material    wallMaterial;   // material name of this wall
    }

    cylinder
    {
        type        foamPatchWall;   // type of the wall

        patch       cylinderWall;       
        
        material    wallMaterial;   // material name of this wall
    }
}
```
 
The most important setup file for CFD-DEM simulation is `constant/couplingProperties`. It contains the parameters for coupling between CFD and DEM, such as drag force closure model, porosity model and etc. It contains two main sub-dictionaries: `unresolved` and `particleMapping`. The `unresolved` dictionary contains the parameters for unresolved coupling, while the `particleMapping` dictionary contains the parameters for particle mapping onto the CFD mesh and MPI parallelization of simulation.

`unresolved` sub-dictionary contains these main parts:

- `distributionMethod`: This parameter defines the method of distributing particle properties (like volume, drag force) across the cells. The options are:

  - `PCM`: Particle Centroid Method - no distribution, direct assignment to cell containing particle center.
  - `Gaussian`: Distributes property over surrounding cells based on a Gaussian distribution using a specified `standardDeviation` value (distribution width).
  - `GaussianIntegral`: Similar to `Gaussian`, but uses the integral of the Gaussian function for distribution.
  - `adaptiveGaussian`: Similar to Gaussian method, but adapts the distribution based on local cell size and particle size. This is the most flexible and accurate method for distributing particle properties across cells.
  - `diffusion`: Uses Laplacian diffusion smoothing to distribute particle properties across cells.
  - `subDivision29`: Divides particle sphere into 29 equal volumetric segments for high-accuracy porosity calculations.
  - `subDivision9`: Divides particle sphere into 9 equal volumetric segments, offering good balance between accuracy and computational cost.

- `porosity`: This part defines the method for calculating porosity. The options are:
  - `distribution`: Uses the selected `distributionMethod` to distribute particle volume over cells and calculate porosity.
  - `subDivision29`: Uses the 29-subdivision method for high-accuracy porosity calculations.
  - `subDivision9`: Uses the 9-subdivision method for porosity calculations.

- `momentumInteraction`: This section contains all momentum coupling related settings:

  - `momentumExchange`: Controls how momentum exchange terms are distributed to fluid cells. Options are:
    - `cell`: No smoothing, assigned to cell containing particle.
    - `distribution`: Uses `distributionMethod` for smoothing momentum terms across cells.

  - `fluidVelocity`: Defines how fluid velocity at particle location is evaluated:
    - `cell`: Uses fluid velocity of the cell containing the particle center.
    - `interpolate`: Interpolates fluid velocity from neighboring cells to particle center.
    - `distribution`: Uses `distributionMethod` for volume-averaged velocity evaluation.

  - `solidVelocity`: Defines how particle velocity is evaluated in coupling calculations:
    - `particle`: Uses exact particle velocity directly.
    - `distribution`: Uses `distributionMethod` to obtain volume-averaged velocity in cell; this average is used for all particles in that cell.

  - `drag`: Defines drag force model and parameters:
    - `model`: Drag closure options are `DiFelice`, `ErgunWenYu`, `Beetstra`, `Rong`, `Cello`.
    - `residualRe`: Minimum Reynolds number threshold to prevent numerical issues.

  - `lift`: Defines lift force model and surface rotation torque:
    - `model`: Lift force options are `none` (default), `Saffmann`, `Loth2008`, `Shi2019`.
    - `surfaceRotationTorque`: Torque model options are `none`, `lowReynolds`, `Loth2008`, `Shi2019`.
    - `residualRe`: Minimum Reynolds number threshold.

For detailed mathematical formulations and more information about the coupling properties, please refer to the main coupling system documentation in `phasicFlowCoupling/couplingSystem/unresolved/README.md`. 