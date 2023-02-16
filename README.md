# PhasicFlow-coupling
PhasicFlow-coupling is a software package for simulating fluid-particle flows. The fluid is assumed as a continuum phase and the particles as discrete bodies. Fluid is modeled using CFD and particles using DEM.

DEM calculations are handled using features of [PhasicFlow](https://github.com/PhasicFlow/phasicFlow). It a parallel DEM package that can be run on multi-core CPUs or GPUs. The equations for the fluid phase are discritized and solved using using OpenFOAM. OpenFOAM is parallelized based on MPI for being executed on multicore CPUs. The fluid-particle coupling uses both parallelization methods, shared-memory and MPI, to leverage the maximum computational resoureces. 

Based on the above configuration, PhasicFlow-coupling can use purly the computational resources of a multi-core CPU or use the those of a computer with both CPU and GPUs. 

## How to build
You need to first install PhasicFlow and OpenFoam on your computer. After that, navigate to the main directory of the code and enter the following command to install the code.

`> ./Allwmake`


