# Pneumatic Conveying in an Elbow Pipe

## 1. Problem Description

This tutorial demonstrates the simulation of a gas flow through a 90-degree elbow pipe using the `unresolvedSpherePFPlus` solver.The particles are inserted into the pipe (diameter is 10 cm) from the horizontal inlet and are transported by the fluid flow (here, it is air). The insertion of particles is done for the first 1 second of the simulation, after which no new particles are added. Air is injected at a superficial velocity of 15 m/s, while the particles are injected at a slower velocity of 4 m/s.

<div align="center">
<b>
<img src="./elbowFlow.gif" alt="Pneumatic conveying in an elbow pipe" />
</b>
<b>

A visualization of a pneumatic conveying in an elbow pipe.
</b></div>


## 2. Running the Simulation

A script named `runThisCase` is provided to automate the process. Below is a step-by-step explanation of the workflow.

### 2.1. Mesh Generation for CFD

The mesh is provided in Ansys Fluent format (`elbow.msh`). It must be converted to OpenFOAM format.

```bash
fluent3DMeshToFoam elbow.msh
```

### 2.2. Particle Phase Initialization

Initialize the discrete particles based on the configuration in `settings/particlesDict`. This step generates the initial particle positions and their fields.

```bash
particlesPhasicFlow
```

### 2.3. Geometry Creation for DEM

Initialize the boundary geometries for the particle solver as defined in `settings/geometryDict`. The surface of the pipe is constructed from the mesh directly. Here, the type should be set to `foamPatchWall` to use the mesh boundary as the wall for the DEM part. `patch` should match the name of the boundary in the mesh (in this case, `shell`), and `material` should be set to the material defined in `caseSetup/interaction` (here, `wallMat`).

```C++
// file: settings/geometryDict

// some lines are missing here
surfaces
{

    shellWall
    {
        type        foamPatchWall;  
        patch       shell;
        material    wallMat;       
    }

}
```

run the following command to create the geometry for DEM part:

```bash
geometryPhasicFlow
```

### 4. Fluid Field Initialization

Copy the initial fluid field files (`U`, `p`, etc.) from the `Fluid` directory to the `0` directory (start time).

```bash
cp -rv ./Fluid/* 0
```

### 5. Running the Solver

Execute the main coupled solver.

```bash
unresolvedSpherePFPlus
```

If you have enough computation resoures, you can run the solver in parallel mode using mpirun. Suppose you want to use 4 processors, the command whould be:

```bash
dcomposePar 
mpirun -np 4 unresolvedSpherePFPlus -parallel
```

Note that in number of processors should match the number of subdomains defined in `system/decomposeParDict`.

### 6. Post-processing

Convert the particle data (lagrangian) to VTK format for easy visualization in Paraview. We export specific fields like velocity, diameter, and ID.

```bash
pFlowToVTK --binary --fields velocity diameter id
```

Create a dummy file to open the OpenFOAM case in Paraview.

```bash
touch foam.foam
```

## Visualization

To view the results:

1. Open **Paraview**.
2. Open `foam.foam` to visualize the Fluid fields (Velocity, Pressure).
3. Open the generated VTK files (in `VTK/` folder usually) to visualize the Particles.
