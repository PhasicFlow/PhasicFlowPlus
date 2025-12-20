# sandBridge

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: unresolvedSpherePFPlus

## 0. Problem Definition

In this tutorial, we will simulate a liquid-solid flow to determine when a bridge is formed in a conduit with narrow gap using the unresolved solver `unresolvedSpherePFPlus`. The geometry has dimensions of 0.077 × 0.091 × 0.007 m³, and spherical particles with a diameter of 0.0023 m and a density of 2500 kg/m³ are simulated. In this simulation, both particles and fluid flow simultaneously. Water is uniformly injected from the upper part of the geometry with a superficial velocity of 0.02 m/s. The simulation runs for a total of 7 seconds.



<div align="center">
<b>
<img src="./sandBridge.gif" alt="sandBridge" style="width: 720px;"/>
</b>
<b>

A visualization of a fluid-particle flow with the fluid and particle field colored based on velocity.
</b></div>

***

## 1. Performing the Simulation Using the Allrun Script

The `Allrun` script is designed to automate the simulation process for the liquid-solid flow for sand bridge formation using the `unresolvedSpherePFPlus` solver. It manages all essential steps, including mesh generation, CFD-DEM coupling simulation, and result conversion. To execute the simulation, follow these steps:

### Step 1: Execute the `Allrun` Script

1. Navigate to the `sandBridge` directory.
2. Run the `Allrun` script by executing the following command:

   ```sh
   ./Allrun
   ```

   This script automates the entire simulation workflow, including mesh generation, CFD-DEM coupling, and result conversion.

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

First part about particle insertion in `caseSetup` folder. This configuration enables particle insertion `(active yes)` in the defined box region called `topRegion`. Particles are inserted from simulation time 0 to 3.5 seconds at a rate of 1000 particles per second, with an insertion interval of 0.049 seconds. The insertion region is a thin box defined by the coordinates min (0.005 0.05 0.002) and max (0.075 0.052 0.005). Inserted particles are given an initial velocity of (0 -0.05 0), meaning they move downward. The particle type used for insertion is `sphere1`, as specified in the mixture section. This setup allows regular particle injection into a narrow area with controlled timing and velocity.

```C++
// This activates particle insertion. If set to no, the particle insertion must be done in paritcleDict
active  yes;

/*
one region is considered for inserting particles. 
*/
topRegion
{
    // Controls insertion based on simulation time
    timeControl  simulationTime;

	// insertion rate (particles/s)
	rate 	 1000;

	// Start time of Particles insertion (s)
	startTime 	  0;

	// End time of Particles insertion (s)	
	endTime   	  7;

	// Time Interval between each insertion (s)
	insertionInterval       0.049;
	
    // type of insertion region and other option is avalible
    regionType   box;
	
    // Coordinates of BoxRegion (m,m,m)
	boxInfo 
	{
		min ( 0.005 0.05 0.002); //  (m,m,m)
		max ( 0.075 0.052 0.005); // (m,m,m)
	}

    setFields
    {
        // initial velocity of inserted particles
        velocity    realx3 (0 -0.05 0); 
    }
   
    mixture
    {
        // first layer of inserted particles and use for binary insertion
        sphere1 1;
    }
}
```
The most important setup file for CFD-DEM simulation is `constant/couplingProperties`. It contains parameters for coupling between CFD and DEM, with two main sub-dictionaries: `unresolved` and `particleMapping`. 

In this tutorial, `distributionMethod` is set to `diffusion` with `standardDeviation=0.0075` (about 3× particle diameter) and `nSteps=10` for smoothing particle properties across cells. The drag model is `DiFelice`, and momentum exchange is distributed across cells using the diffusion method.

```C++
// constant/couplingProperties file 
unresolved
{
    // Distribution method for mapping particle data across fluid cells
    distributionMethod      diffusion;
    
    // Distribution method settings 
    diffusionInfo
    {
        nSteps              10;
        standardDeviation   0.0075;
        
        log                 0;
    }

    // Porosity calculation method 
    porosity
    {
        method      distribution;
        alphaMin    0.2;
    }

    // Momentum coupling settings
    momentumInteraction
    {
        momentumExchange distribution; 
        fluidVelocity    cell;
        solidVelocity    particle;

        drag
        {
            model       DiFelice; 
            residualRe  1.0e-6;
        }

        lift
        {   
            model                   none;
            surfaceRotationTorque   none; 
            residualRe              1.0e-6;
        }

        virtualMass
        {
            // Not implemented yet
        }
    }

}

particleMapping
{
    domainExpansionRatio    1;
    domainUpdateInterval    0.0001;
    decompositionMode       facePlanes;
}
```

For detailed information about coupling parameters, refer to the main coupling system documentation.
