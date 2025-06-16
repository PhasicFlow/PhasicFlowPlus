# Fluidized bed

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: unresolvedSpherePFPlus

## 0. Problem Definition

In this tutorial, we will simulate a gas-solid fluidized bed using the unresolved solver `unresolvedSpherePFPlus`. The fluidized bed has dimensions of 0.15x0.7x0.04 m³ and contains 100,000 spherical particles with a diameter of 0.0018 m and a density of 1000 kg/m³. Initially, the particles are at rest. Gas is uniformly injected from the bottom of the bed at a superficial velocity of 1.3 m/s. The simulation runs for a total of 10 seconds, with the first second dedicated to the initial packing of particles (pure DEM simulation) and the remaining 9 seconds focused on the fluidized bed simulation.

<div align="center">
<b>
<img src="./fluididze-bed-cfd-dem.jpeg" alt="Fluidized bed" style="width: 400px;"/>
</b>
<b>

A visualization of a gas-solid fluidized bed with the gas field colored based on velocity.
</b></div>

***

## 1. Performing the Simulation Using the Allrun Script

The `Allrun` script is designed to automate the simulation process for the gas-solid fluidized bed using the `unresolvedSpherePFPlus` solver. It manages all essential steps, including mesh generation, DEM simulation, CFD-DEM coupling simulation, and result conversion. The first second of the simulation is dedicated to particle settling, which is a pure DEM simulation using the `sphereGranFlow` solver. Following this phase, the CFD-DEM simulation is performed for the remaining 9 seconds. To execute the simulation, follow these steps:

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

To learn about how to set up a DEM simulation, please refer to the [tutorial page](https://github.com/PhasicFlow/phasicFlow/wiki/Tutorials) of PhasicFlow and other online documents along side this package. Also, you can refer to OpenFOAM tutorials to learn about how to set up a CFD simulation. The solver we are using here, essentially is a combination of DEM and CFD component, with some additional parameters that are essential for unresolved coupling. Here, we only describe the simulation setup files that are specific to the coupling part.

The most important setup file for CFD-DEM simulation is `constant/couplingProperties`. It contains the parameters for coupling between CFD and DEM, such as drag force closure model, porosity model and etc. It contains two main sub-dictionaries: `unresolved` and `particleMapping`. The `unresolved` dictionary contains the parameters for unresolved coupling, while the `particleMapping` dictionary contains the parameters for particle onto the CFD mesh and MPI parallelization of simulation. The parameters in `particleMapping` sub-dictionary are left unchanged, since they are the best settings and we rarely want to change them, although they are here to enable users to manage some special cases.

`unresolved` sub-dictionary contains these parts:

- `cellDistribution`: This part defines the method of distributing particle properties (like volume, drag force) across the cells. The options are `self`, `Gaussian`, `GassianIntegral`, and `adaptiveGaussian`:

  - `self`: distributes property over the cell on which the center of particles is located.
  - `Gaussian`: distributes property over the surrounding cells based on a Gaussian distribution, using a specified `standardDeviation` value (distribution width).
  - `GaussianIntegral`: similar to `Gaussian`, but it uses the integral of the Gaussian function for distribution.
  - `adaptiveGaussian`: method is used, which adapts the distribution based on the local cell size and particle size. This is the most flexible and accurate statistical method for distributing particle properties across the cells.

- `porosity`: This part defines the method for calculating porosity. The options are `PIC`, `subDivision29Mod`, `subDivision9`, `diffusion`, and `cellDistribution`:

  - `PIC`: Particle-In-Cell method.
  - `subDivision29`, `subDivision29Mod` and `subDivision9`: These methods are used for calculating porosity based on the subdivision of particles into equal volumes and mapping these sub-volumes onto the surrounding cells for calculating porosity.
  - `diffusion`: This method uses a diffusion model to calculate porosity.
  - `cellDistribution`: This method uses the cell distribution function defined in the `cellDistribution` sub-dictionary to distribute particle volume over cells and finally calculate porosity.

- `drag`: This part defines the drag force closure model. The options are `DiFelice`, `ErgunWenYu`, and `Rong`. 
  - `fluidVelocity`: This parameter defines how the fluid velocity at the center point of the particle is calculated:
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
