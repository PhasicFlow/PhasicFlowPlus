# Fluidized bed

## Problem definition 

We are going to simulate a gas-solid fluidized bed using the unresolved solver `unresolvedSpherePFPlus`. The dimensions of the fluidized bed is 0.15x0.7x0.04 m3. The bed is filled with 100,000 spherical particles of 0.0018 m diameter with density 1000 kg/m3, which are initially at rest. The gas is uniformly injected from the bottom of the bed with a superficial velocity of 1.3 m/s. The simulation will run for 10 seconds. The first second of the simulation is used for initial packing of particles (pure DEM simulation) and the following 9 seconds are used for the fluidized bed simulation.

<div align="center">
<b>
<img src="./fluididze-bed-cfd-dem.jpeg" alt="Fluidized bed" style="width: 400px;"/>
</b>
<b>

A view of gas-solid fluidized bed with gas field colored based on the velocity
</b></div>

***

## 1. How to Perform the Simulation using Allrun Script

`Allrun` is a script designed to automate the simulation process for the gas-solid fluidized bed using the `unresolvedSpherePFPlus` solver. It handles all necessary steps, including mesh generation, DEM simulation, CFD-DEM coupling simulation, and result conversion. , 
To simulate the gas-solid fluidized bed, follow these steps:

### Step 1: Run the `Allrun` Script

1. Navigate to the `fluidizedbed` directory.
2. Execute the `Allrun` script:
   ```sh
   ./Allrun
   ```

   This script automates the simulation process, including mesh generation, DEM simulation, CFD-DEM coupling, and result conversion.

### Step 2: Understand the Sub-Folders

The first second of simulation is dedicated to the settling of particles. It is a pure DEM simulation with solver `sphereGranFlow`. After this phase, the CFD-DEM simulation is pefromed up to 10 s. The folder structure of simulation consistes of two main parts, folders that contain files related to DEM parameters, and folder that contain files related to CFD and coupling parameters:

- **DEM related folders**:
  - **`settings/`**: Contains configuration files for the DEM simulation.
  - **`caseSetup/`**: Includes files for setting up the simulation/physical parameters for particles.
- **CFD-DEM related folders**:
  - **`FluidField/`**: Holds the initial field data (`alpha`, `p`, `U`) used for the CFD-DEM simulation.
  - **`constant/`**: Contains constant properties for fluid and parameters for coupling (CFD-DEM).
  - **`system/`**: Contains files for setting up CFD simulation parameters.

### Step 3: View Results

After the simulation completes, results are converted to VTK format for visualization. You can find the VTK files in the `./VTK` folder. Execute the following command to visualize the reusults:
```sh
paraview foam.foam &
```
use the `foam.foam` file to open the CFD results in ParaView. And to visualize DEM results, open `./VTK/particles.vtk.series`

## 2. Description of setup files

The most important setup file for CFD-DEM simulation is constant/couplingProperties. It contains the parameters for coupling between CFD and DEM, such as drag force closure model, porosity model and etc. It contains two main sub-dictionaries: `unresolved` and `particleMapping`. The `unresolved` dictionary contains the parameters for unresolved coupling, while the `particleMapping` dictionary contains the parameters for particle onto the CFD mesh and MPI parallelization of simulation. The parameters `particleMapping` sub-dictionary are left unchanged, since they are the best settings and we rarely want to change them, although they are here enable user to manage some special cases.

`unresolved` contains some parts:

- `cellDistribution`: This part defines the method of distributing particle properties (like volume, drag force) across the cells. The options are `self`, `Gaussian`, `GassianIntegral`, and `adaptiveGaussian`:

  - `self`: distributes property over the cell on which the center of particles is located. 
  - `Gaussian`: distributes property over the surrounding cells based on a Gaussian distribution, using a specified `standardDeviation` value (distribution width).
  - `GaussianIntegral`: similar to `Gaussian`, but it uses the integral of the Gaussian function for distribution.
  - `adaptiveGaussian`: method is used, which adapts the distribution based on the local cell size and particle size. This is the most felxible and accurate statistical method for distributing particle properties across the cells.

- `porosity`: This part defines the method for calculating porosity. The options are `PIC`, `subDivision29Mod`, `subDivision9`, `diffusion`, and `cellDistribution`:

  - `PIC`: Particle-In-Cell method.
  - `subDivision29`, `subDivision29Mod` and `subDivision9`: These methods are used for calculating porosity based on the subdivision of particles into equal volumes and mapping these sub-volumes onto the surrounding cells for calculating porosity.
  - `diffusion`: This method uses a diffusion model to calculate porosity.
  - `cellDistribution`: This method uses cell distribution function defined in `cellDistribution` sub-dictionary to distribute particle volume over cells and finally calculating porosity.

- `drag`: This part defines the drag force closure model. The options are `DiFelice`, `ErgunWenYu`, and `Rong`. 
  - `fluidVelocity`: This parameter defines how the fluid velocity at center point of particle is calculated:
    - `cell`: Uses the fluid velocity of the cell that contains the particle center.
    - `interpolation`: Uses averaged fluid velocity based on cell velocities around the particle.
  - `cellDistribution`: This parameter defines whether the calculated drag force is distributed on cells or not. The options are `on` and `off`.


```C++
// constant/couplingProperties file 
unresolved
{

    cellDistribution
    {
        // type of cell distribution method (if required) 
        //    self: no distribution (cell itself)
        //    Gaussian: distribute values on sorounding cells based on a neighbor length
        //    adaptiveGaussian: similar to Gaussian, but it adapts the distribution 
        //    GaussianIntegral 
        type                adaptiveGaussian; 
    }

    porosity
    {
    	// Options are PIC, subDivision29Mod, subDivision9, diffusion, cellDistribution
        method      cellDistribution;

        // minimum alpha allowed 
        alphaMin    0.2;
    }

    drag
    {
        // Drag force closure, other options are ErgunWenYu, Rong
        type                DiFelice; 

        // type of fluid velocity used in the drag force calculations
        //   cell: uses fluid velocity of the cell that contains the particle center 
        //   interpolation: uses averaged fluid velocity based on cell velocities around particle
        fluidVelocity       cell;  

        // weather to distribute calculated drag force on cells
        //   off: set the calculated drag force on cell itself
        //   on: distributed the calculated drag force on cells (using cellDistribution method)
        cellDistribution    off; 

        residualRe          10e-6;
    }

}

particleMapping
{
	// based on the maximum particle diameter in the simulation.
    domainExpansionRatio    1;

    domainUpdateInterval    0.01;

    decompositionMode       facePlanes;
}
```
