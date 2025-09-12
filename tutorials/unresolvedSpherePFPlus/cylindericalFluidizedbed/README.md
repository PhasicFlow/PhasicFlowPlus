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
 
The most important setup file for CFD-DEM simulation is `constant/couplingProperties`. To learn more about this file and how to set it up, you are refered to other tutorials, in which this has been completely explained. 