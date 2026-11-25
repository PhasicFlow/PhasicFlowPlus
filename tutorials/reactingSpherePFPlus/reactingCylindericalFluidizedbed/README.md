# Cylinderical Fluidized Bed (Reacting)

    - Compatibility: OpenFOAM v2406 and v2412 and PhasicFlow-v-1.0.
    - Solver: reactingSpherePFPlus

## 0. Problem Definition

In this tutorial, we simulate the fuel reactor of chemical looping
combustion: a batch fluidized bed of oxygen carrier that is reduced by
carbon monoxide,

3 Fe2O3 + CO → 2 Fe3O4 + CO2

using the reacting solver `reactingSpherePFPlus`. The column is 0.05 m in
diameter and 0.3 m high (the
[thermal tutorial](../../unresolvedHeatSpherePFPlus/heatCylindericalFluidizedbed)
has 0.1 m and 0.2 m): narrower, so that the run stays light, and taller, to
leave room for the bubbling bed. It holds 12.5K porous particles of 2 mm,
which settle into a bed about 4.6 cm high. They are the carrier Fe45Al-FG of
Abad et al. (2007): 45 wt% active Fe2O3, the rest Al2O3 and inactive Fe2O3.
The bed, the gas and the furnace wall are at 1173.15 K (900 °C), at 1 atm.

The run has three stages:

- 0 to 0.25 s: the particles settle, DEM only (`multiSpeciesGranFlow`), in
  N2, so nothing reacts.
- 0.25 to 1 s: N2 enters from the bottom at 3.0 m/s (about 2.2 times the
  minimum fluidization velocity) and fluidizes the bed.
- 1 to 5 s: the inlet gas switches to 50 vol% CO in N2. The CO reduces the
  particles, and CO2 leaves the bed with the unconverted CO.

The run ends partway through the reduction: the CO that the whole bed needs
enters in about 5 s, but a single particle takes about 35 s at the inlet CO,
since the gas has to diffuse into its pores. Raise `endTime` in
`system/controlDict` for a longer run.

***

## 1. Performing the Simulation Using the Allrun Script

The `Allrun` script runs all the steps:

1. `blockMesh` builds the mesh.
2. `particlesPhasicFlow`, `geometryPhasicFlow`, `correctionTablePhasicFlow`
   and `multiSpeciesGranFlow` place the particles, build the walls, build
   the correction table of the grain model and run the settling stage
   (0 to 0.25 s).
3. The fluid fields (`alpha`, `p`, `U`, `T`, `Y.CO`, `Y.CO2`, `Y.N2`) are
   copied from `FluidField/` to the time folder 0.25.
4. `reactingSpherePFPlus` runs the CFD-DEM stage (0.25 to 5 s).
5. `pFlowToVTK` converts the particle results for ParaView.

```sh
./Allrun
```

## 2. What's New Compared With the Thermal Tutorial

Please go through the thermal tutorial first. The setup files below are
the ones that change; the sources of all values are given next to them.

### Gas species

`constant/transportProperties` lists the gas species and their properties.
The solver computes the mixture density from the ideal-gas law at `pRef`,
and the mixture heat capacity and conductivity from the mass fractions:

```C++
gasSpecies      (CO      CO2     N2);
gasMw           (28.010  44.010  28.013);   // [kg/kmol]
gasCp           (1208    1265    1194);     // [J/(kg.K)]
gasKappa        (0.0740  0.0785  0.0742);   // [W/(m.K)]
inertSpecies    N2;                         // not transported: 1 - the others
Sc              [0 0 0 0 0 0 0] 0.666;      // D = nu/Sc
pRef            [1 -1 -2 0 0 0 0] 101325;
```

Each species has a field `Y.<name>` (mass fraction) in `FluidField/`. The
inlet value of `Y.CO` and `Y.N2` is a table in time: N2 until 1 s, then 50
vol% CO. `system/fvSchemes` has a divergence scheme for each transported
species, and `system/fvSolution` solves them with the `Y.*` entries.

### Time-step limit of the reactions

`system/fvSolution` has a `CLC` dictionary. The time step is lowered when
the gas mass of a cell, or the amount of one consumed gas, would change by
more than these fractions in one step (`adjustTimeStep yes` in
`system/controlDict`):

```C++
CLC
{
    maxReactionNo           0.15;   // gas mass of a cell
    maxSpeciesReactionNo    0.15;   // each consumed gas
}
```

### Reaction coupling

`constant/couplingProperties` gains a `reactionInteraction` block under
`unresolved`:

```C++
reactionInteraction
{
    reactionExchange    distribution;   // cell | distribution

    reactionTransferModel
    {
        model           sphere;
        productGas      (CO CO2 N2);     // same order as gasSpecies
    }

    massTransferModel
    {
        model           Gunn;            // RanzMarshall | Gunn | none
    }
}
```

- `reactionExchange` maps the gas sources of the particles onto the mesh,
  like `heatSourceExchange` for heat.
- `productGas` names the gas species that the DEM reactions exchange; it
  must match `gasSpecies` and the `productGas` of the DEM.
- `massTransferModel` computes the gas-film coefficient k_c of each particle
  from one of the Nusselt closures of `heatTransferModel`, with the Schmidt
  number in place of the Prandtl number, and sends it to the DEM. Without
  this block the DEM keeps the `k_c` of its reaction dictionary.

`heatTransferModel` uses the `Gunn` closure here, which depends on the fluid
volume fraction and suits a dense bed.

### Reactions of the particles

`caseSetup/thermoPhysicalInteraction` holds the reactions, as in the
`packedBedOxidation` tutorial of phasicFlow: the solid and gas species, the
kinetics of Abad et al. (2007), the grain model of the porous particle with
its correction table (`correctionTable Yes`, built by
`correctionTablePhasicFlow`), and `defaultGasConcentrations`, the gas that
the particles see before the first CFD exchange. `settings/particlesDict`
sets the `solidFrac_<species>` fields of the particles.

### Heat transfer

All particle heat paths are on: particle-particle radiation, contact
conduction and conduction through the gas between particles, and the three
particle-wall paths. The inlet and the cylinder are at 1173.15 K (key
`temperature` in `settings/geometryDict`); the outlet has no temperature and
stays adiabatic. The wall material is Kanthal APMT, the FeCrAl alloy of a lab
CLC reactor. For the gas, the wall patch has `fixedValue` 1173.15 K in
`FluidField/T`.

## 3. Results

- `postProcessing/outletGas/`: flux-weighted mean of `Y.CO`, `Y.CO2`, `Y.N2`
  and `T` over the outlet, every 0.01 s. The CO2 that leaves the bed is the
  measure of the reduction.
- `postProcessing/bedGas/`: volume-mean gas temperature of the column.
- `postprocessData/`: history of 10 particles (ids 0, 1250, ..., 11250)
  every 0.05 s: `temperature`, `gasConc_CO` (the CO concentration the
  particle sees, kmol/m³), `solidFrac_Fe2O3` and `massTransferCoeff` (k_c).
  `gasConc_<gas>` and `massTransferCoeff` live in memory only, so probing
  them here is the way to record them.
- `wallHeatFlow.dat`: the heat that each wall surface gives to the
  particles through each path.

To visualize the results:

```sh
paraview foam.foam &
```

For the particles, open `./VTK/particles.vtk.series`.

## References

- A. Abad, J. Adánez, F. García-Labiano, L.F. de Diego, P. Gayán, J. Celaya,
  Mapping of the range of operational conditions for Cu-, Fe-, and Ni-based
  oxygen carriers in chemical-looping combustion, *Chem. Eng. Sci.* 62
  (2007) 533-549: carrier, Tables 2 and 3.
- F. García-Labiano, L.F. de Diego, J. Adánez, A. Abad, P. Gayán, Temperature
  variations in the oxygen carrier particles during their reduction and
  oxidation in a chemical-looping combustion system, *Chem. Eng. Sci.* 60
  (2005) 851-862: effective diffusivity, particle conductivity and
  emissivity.
- D.J. Gunn, Transfer of heat or mass to particles in fixed and fluidised
  beds, *Int. J. Heat Mass Transfer* 21 (1978) 467-476.
- C.Y. Wen, Y.H. Yu, A generalized method for predicting the minimum
  fluidization velocity, *AIChE J.* 12 (1966) 610-612.

For the coupling formulations, see
`phasicFlowCoupling/couplingSystem/unresolved/README.md`.
