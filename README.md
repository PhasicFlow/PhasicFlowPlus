# PhasicFlowPlus
PhasicFlowPlus is a software package for simulating fluid-particle flows. It is a combination of computational fluid dynamics (CFD) and discrete element method (DEM). The fluid is assumed as a continuum phase and the particles as discrete bodies.

DEM calculations are handled using features of [PhasicFlow](https://github.com/PhasicFlow/phasicFlow). It is a parallel DEM package that can be run on multi-core CPUs or GPUs. The equations for the fluid phase are discritized based on finite volume (FV) and solved using using OpenFOAM. OpenFOAM is parallelized based on MPI for being executed on multicore CPUs. The fluid-particle coupling uses both parallelization methods, shared-memory and MPI, to leverage the maximum computational resoureces. 

Based on the above configuration, PhasicFlowPlus can use the computational resources of a multi-core CPU or use the those of a computer with both CPU and GPUs. 

## How to build
You need to first install PhasicFlow and OpenFoam (For now, it is only tested with [OpenFOAM-v9](https://openfoam.org/download/9-source/)) on your computer. After that, clone PhasicFlowPlus on your computer and navigate to the root directory of the code and enter the following command to install the code.

`> ./Allwmake`


