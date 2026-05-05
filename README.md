# laserHeatFoam

`laserHeatFoam` is an OpenFOAM solver for Laser Powder Bed Fusion
(L-PBF) thermal simulations.

It solves transient heat conduction with:

- a moving volumetric laser source,
- temperature-dependent material properties,
- latent heat via apparent heat capacity,
- surface loss from radiation, convection, and evaporation,
- layer-by-layer activation with `fvMeshSubset`,
- optional mesh refinement with `dynamicRefineFvMesh`.

The merged PRs #2 and #3 added the multi-layer workflow, the
`tutorials/multiLayer_basic` case, and the AMR path used by that
tutorial.

## Key Files

- `applications/solvers/laserHeatFoam/laserHeatFoam.C`
- `applications/solvers/laserHeatFoam/createFields.H`
- `applications/solvers/laserHeatFoam/multiLayer.H`
- `applications/solvers/laserHeatFoam/updateRefineIndicator.H`
- `applications/solvers/laserHeatFoam/rebuildActiveCells.H`
- `tutorials/multiLayer_basic`

## Multi-Layer Mode

Set the layer-activation controls in
`constant/multiLayerProperties`.

- `baseplatePatch`: patch that seeds layer 0.
- `nLayers`: number of layers to activate.
- `cellsPerLayer`: cells grown per layer in the build direction.
- `layerDuration`: simulated time assigned to each layer.
- `layerInitialT`: temperature applied when a layer is activated.
- `exposedFacesPatch`: optional patch for the exposed subset faces.

Layer 0 is formed from cells touching `baseplatePatch`.
Later layers are grown by face-neighbour expansion.
The solver keeps the base mesh fields, solves on a sub-mesh, and then
scatters the results back to the base mesh for writing.

## AMR

If `constant/dynamicMeshDict` is present, the solver can refine the base
mesh near the melt-history interface.

- `dynamicFvMesh dynamicRefineFvMesh`
- `field refineIndicator`
- `lowerRefineLevel` and `upperRefineLevel`
- `unrefineLevel`
- `indicatorSmoothing`

The `layerID` field records layer membership.
The `refineIndicator` field drives refinement.

## Case Setup

Required inputs:

- `constant/LaserProperties`
- `constant/timeVsLaserPosition`
- `constant/timeVsLaserPower`
- `constant/thermophysicalProperties`
- `constant/multiLayerProperties`
- `0/T`
- `0/rho`
- `0/cp`
- `0/Q`
- `0/rc`

The multi-layer tutorial is
`tutorials/multiLayer_basic`.

## Build

Use OpenFOAM v2512.
Source the environment before building.

```bash
source ~/bin/load-openfoam v2512
cd applications/boundaryConditions/pbfRadEvapTemperature
wmake libso
cd ../../solvers/laserHeatFoam
wmake
```

## Run

Single case:

```bash
cd tutorials/LPBF_titanium_case_195W_1mps
./Allrun
```

Multi-layer tutorial:

```bash
cd tutorials/multiLayer_basic
./Allrun
```

Parallel tutorial runs use `runMode=parallel`.
Serial runs use `runMode=serial`.

## Output Fields

The solver writes these fields:

- `T`
- `Q`
- `k`
- `rc`
- `fL`
- `rp`
- `rm`
- `rs`
- `rho`
- `cp`
- `cpEff`
- `meltHistory`
- `layerID`
- `refineIndicator`
- `gradT`

## References

- Mohammadkamal and Caiazzo (2025),
  [The role of laser operation mode on thermal and mechanical behavior
  in powder bed fusion](https://doi.org/10.1007/s00170-025-16460-4)
- Proell et al. (2023),
  [A highly efficient computational approach for fast scan-resolved
  simulations of metal additive manufacturing processes on the scale of
  real parts](https://arxiv.org/abs/2302.05164)
- Proell et al. (2020),
  On phase change and latent heat models in metal additive manufacturing
  process simulation.
- Anisimov and Khokhlov (1995),
  Instabilities in laser-matter interaction.

## License

This solver is distributed under the GNU General Public License v3.
