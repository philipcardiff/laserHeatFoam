# laserHeatFoam

An OpenFOAM solver for thermal simulation of **Laser Powder Bed Fusion (L-PBF)** additive manufacturing processes. The solver resolves the moving laser heat source, phase-dependent material properties, latent heat of fusion, physically consistent surface boundary conditions including radiation, conduction and evaporation, and progressive layer-by-layer activation for multi-layer cases.

The physical and numerical models are based on:

- **[Paper 1]** Mohammadkamal & Caiazzo (2025), *"The role of laser operation mode on thermal and mechanical behavior in powder bed fusion: a numerical study"*, The International Journal of Advanced Manufacturing Technology, 140:3779–3796. https://doi.org/10.1007/s00170-025-16460-4
- **[Paper 2]** Proell, Munch, Kronbichler, Wall & Meier (2023), *"A highly efficient computational approach for fast scan-resolved simulations of metal additive manufacturing processes on the scale of real parts"*, arXiv:2302.05164v3.

---

## Table of Contents

1. [Governing Equation](#1-governing-equation)
2. [Volumetric Heat Source Models](#2-volumetric-heat-source-models)
3. [Phase Change and Latent Heat Model](#3-phase-change-and-latent-heat-model)
4. [Phase-Dependent Thermal Conductivity](#4-phase-dependent-thermal-conductivity)
5. [Temperature-Dependent Material Properties](#5-temperature-dependent-material-properties)
6. [Boundary Conditions](#6-boundary-conditions)
7. [Diffusion Number and Time Step Control](#7-diffusion-number-and-time-step-control)
8. [Solver Algorithm](#8-solver-algorithm)
9. [Multi-Layer Activation and AMR](#9-multi-layer-activation-and-amr)
10. [Case Setup](#10-case-setup)
11. [Building and Running](#11-building-and-running)
12. [Output Fields](#12-output-fields)
13. [References](#13-references)

---

## 1. Governing Equation

The solver solves the **transient heat conduction equation** with a volumetric source term:

$$\rho c_{p}^{\text{eff}} \frac{\partial T}{\partial t} = \nabla \cdot \left( k \nabla T \right) + Q$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $T$ | Temperature | K |
| $\rho$ | Density | kg/m³ |
| $c_p^{\text{eff}}$ | Effective specific heat (includes latent heat, see §3) | J/(kg·K) |
| $k$ | Thermal conductivity (phase-weighted, see §4) | W/(m·K) |
| $Q$ | Volumetric laser heat source (see §2) | W/m³ |

This is implemented in OpenFOAM as an `fvScalarMatrix`:

```cpp
fvScalarMatrix TEqn
(
    rho * cpEff * fvm::ddt(T) - fvm::laplacian(k, T) == Q
);
```

The spatial discretization uses a standard finite-volume approach on an unstructured mesh. Time integration uses the implicit Euler scheme (backward differencing), which is unconditionally stable and suitable for the diffusion-dominated problem.

---

## 2. Volumetric Heat Source Models

The laser is represented as a **moving volumetric heat source** with a Gaussian radial profile. The laser position and power are prescribed as time series in `constant/timeVsLaserPosition` and `constant/timeVsLaserPower`, allowing arbitrary scan paths and continuous or pulsed-wave operation.

Two depth profiles are available via the `laserModel` keyword in `constant/LaserProperties`.

### 2.1 Cylindrical Model (`laserModel cylindrical`)

Based on Proell et al. (2023), eq. (6), and Mohammadkamal & Caiazzo (2025), eq. (4):

$$Q(\hat{x}, \hat{y}, \hat{z}) = \frac{2 \, A \, P_{\text{eff}}}{\pi \, R^2 \, d} \exp\left( -\frac{2(\hat{x}^2 + \hat{y}^2)}{R^2} \right), \quad \text{if } z_{\text{layer}} - d < z < z_{\text{layer}}$$

The source is **uniform in depth** within a cylinder of depth $d$ below the laser surface position $z_{\text{layer}}$, and zero elsewhere.

### 2.2 Exponential Decay Model (`laserModel exponential`)

An alternative model with **Beer–Lambert exponential attenuation** in the $z$-direction:

$$Q(\hat{x}, \hat{y}, z) = \frac{2 \, A \, P_{\text{eff}}}{\pi \, R^2 \, d} \exp\left( -\frac{2(\hat{x}^2 + \hat{y}^2)}{R^2} \right) \exp\left( -\frac{z_{\text{layer}} - z}{d} \right), \quad z < z_{\text{layer}}$$

This is better suited when the laser penetration depth follows an exponential decay law (e.g., for semi-transparent or highly scattering powder beds).

### Heat Source Parameters

| Parameter | Symbol | Dictionary key | Description |
|-----------|--------|---------------|-------------|
| Laser power | $P$ | `timeVsLaserPower` | Time-varying power [W] |
| Absorptivity | $A$ | `absorptivity` | Fraction of power absorbed |
| Beam radius | $R$ | `laserRadius` | 1/e² Gaussian radius [m] |
| Absorption depth | $d$ | `absorptionDepth` | Penetration depth [m] |
| Laser position | — | `timeVsLaserPosition` | Moving center (x, y, z) [m] |

The Gaussian profile in the x-y plane has a standard deviation $\sigma = R/2$, meaning $R$ is the effective beam radius at which the intensity falls to $e^{-2}$ of its peak.

**Pulsed-wave (PW) operation** is naturally supported by specifying a time series for laser power that alternates between the peak power and zero, following eq. (8) of Paper 1:

$$g(t) = \begin{cases} 1, & \text{during pulse-on time} \\ 0, & \text{during pulse-off time} \end{cases}$$

so that $Q_{\text{PW}} = Q_{\text{CW}} \cdot g(t)$. The duty cycle $DC = f \cdot T_p$ (frequency × pulse width) controls the average heat input.

---

## 3. Phase Change and Latent Heat Model

The solver implements the **apparent heat capacity** (enthalpy) method to account for the latent heat of fusion at the solid–liquid interface, as described in Proell et al. (2020, 2023).

### 3.1 Liquid Fraction

The liquid fraction $f_L$ ramps linearly in the mushy zone between the solidus temperature $T_s$ and the liquidus temperature $T_l$:

$$f_L(T) = \begin{cases} 0, & T < T_s \\ \dfrac{T - T_s}{T_l - T_s}, & T_s \leq T \leq T_l \\ 1, & T > T_l \end{cases}$$

### 3.2 Effective Heat Capacity

The latent heat $L$ [J/kg] is absorbed over the mushy zone by augmenting the specific heat:

$$c_p^{\text{eff}} = c_p + L \frac{\partial f_L}{\partial T} = c_p + \frac{L}{T_l - T_s} \cdot \mathbf{1}_{[T_s, T_l]}(T)$$

where $\mathbf{1}_{[T_s, T_l]}$ is an indicator function that is 1 inside the mushy zone and 0 otherwise. This ensures that the correct amount of latent heat is released/absorbed as the material passes through the phase-change interval.

### 3.3 Consolidated Fraction (Irreversible Powder-to-Melt Transition)

Following Proell et al. (2023), eq. (3), the **consolidated fraction** $r_c$ tracks the maximum liquid fraction ever reached by each material point:

$$r_c(t) = \max_{\tilde{t} < t} f_L(T(\tilde{t}))$$

This is an **irreversible** quantity: once powder melts, $r_c$ can only increase. It distinguishes between:
- **Powder** that has never melted ($r_c = 0$, initialized for powder cells)
- **Consolidated material** that has melted at least once ($r_c > 0$)
- **Substrate** that starts consolidated ($r_c = 1$, initialized for solid cells)

In code:
```cpp
rc = max(rc, fL);   // monotonically non-decreasing
```

---

## 4. Phase-Dependent Thermal Conductivity

Following Proell et al. (2023), eq. (5), the thermal conductivity is computed as a weighted average over the three material phases:

$$k(T, r_c) = r_p(r_c)k_p + r_m(T)k_m + r_s(T,r_c)k_s$$

where the phase fractions (Proell 2023, eqs. 4) are:

| Fraction | Formula | Description |
|----------|---------|-------------|
| $r_p = 1 - r_c$ | Powder fraction | Decreases irreversibly once material melts |
| $r_m = f_L$ | Melt fraction | Liquid portion at current temperature |
| $r_s = r_c - f_L$ | Solid fraction | Consolidated but solidified material |

And the single-phase conductivities are:

| Symbol | Description |
|--------|-------------|
| $k_p$ | Powder thermal conductivity (low, due to porosity) |
| $k_m$ | Melt thermal conductivity |
| $k_s$ | Solid thermal conductivity (temperature-dependent) |

The very low powder conductivity ($k_p \ll k_s$) effectively acts as a thermal insulation boundary around the part — physically justified because the powder bed has very poor heat conduction compared to the consolidated material.

---

## 5. Temperature-Dependent Material Properties

For **Ti-6Al-4V**, piecewise polynomial fits are used for density, specific heat, and conductivity, based on Mohammadkamal & Caiazzo (2025), Table 1.

### 5.1 Density

$$\rho_{\text{solid}}(T) = \begin{cases} \rho_{0s} + \rho_{1s} \, T, & T < T_l \\ \rho_{0m} + \rho_{1m} \, T, & T \geq T_l \end{cases}$$

For powder cells, the effective density is reduced by the initial powder porosity $\varphi$ (Mohammadkamal eq. 1):

$$\rho_{\text{eff}} = \rho_{\text{solid}} \cdot (1 - \varphi \cdot r_p)$$

so that $r_p = 1$ (pure powder) gives $\rho_{\text{eff}} = (1-\varphi)\,\rho_{\text{solid}}$ and $r_p = 0$ (solid/melt) gives $\rho_{\text{eff}} = \rho_{\text{solid}}$.

### 5.2 Specific Heat

$$c_p(T) = \begin{cases} c_{p,0a} + c_{p,1a} \, T, & T < T_1 \\ c_{p,0b} + c_{p,1b} \, T, & T_1 \leq T < T_l \\ c_{p,L}, & T \geq T_l \end{cases}$$

where $T_1 = 1268\,\text{K}$ is a breakpoint temperature for Ti-6Al-4V.

### 5.3 Solid/Melt Thermal Conductivity

$$k_{s/m}(T) = \begin{cases} k_{0a} + k_{1a} \, T, & T < T_1 \\ k_{0b} + k_{1b} \, T, & T_1 \leq T < T_l \\ k_{0m} + k_{1m} \, T, & T \geq T_l \end{cases}$$

### Default Coefficients for Ti-6Al-4V

| Coefficient | Value | Units |
|-------------|-------|-------|
| $\rho_{0s}$ | 4557.892 | kg/m³ |
| $\rho_{1s}$ | −0.154 | kg/m³/K |
| $\rho_{0m}$ | 5227.640 | kg/m³ |
| $\rho_{1m}$ | −0.680 | kg/m³/K |
| $c_{p,0a}$ | 483.04 | J/kg/K |
| $c_{p,1a}$ | 0.215 | J/kg/K² |
| $c_{p,0b}$ | 412.70 | J/kg/K |
| $c_{p,1b}$ | 0.18 | J/kg/K² |
| $c_{p,L}$ | 831.0 | J/kg/K |
| $k_{0a}$ | 1.250 | W/m/K |
| $k_{1a}$ | 0.015 | W/m/K² |
| $k_{0b}$ | 3.15 | W/m/K |
| $k_{1b}$ | 0.012 | W/m/K² |
| $k_{0m}$ | −12.752 | W/m/K |
| $k_{1m}$ | 0.024 | W/m/K² |
| $T_s$ | 1877 | K |
| $T_l$ | 1923 | K |
| $T_1$ | 1268 | K |
| $L$ | 286,000 | J/kg |
| $\varphi$ | 0.35 | — |

---

## 6. Boundary Conditions

### 6.1 Top Surface — Radiation, Convection, and Evaporation (`pbfRadEvapTemperature`)

The top free surface of the powder bed / built part is subject to heat loss from three mechanisms, implemented in the custom boundary condition `pbfRadEvapTemperatureFvPatchScalarField`.

The total surface heat flux out of the domain is:

$$q_{\text{tot}} = q_{\text{rad}} + q_{\text{conv}} + q_{\text{evap}}$$

and enters the solver as a Neumann condition:

$$-k \frac{\partial T}{\partial n} \bigg|_{\text{top}} = q_{\text{tot}}$$

#### Radiation (Proell 2023, eq. 11; Mohammadkamal eq. 10)

$$q_{\text{rad}} = \varepsilon \sigma_{\text{SB}} \left( T_w^4 - T_\infty^4 \right)$$

| Symbol | Description |
|--------|-------------|
| $\varepsilon$ | Surface emissivity (Ti-6Al-4V ≈ 0.3–0.4) |
| $\sigma_{\text{SB}}$ | Stefan–Boltzmann constant = 5.67×10⁻⁸ W/m²/K⁴ |
| $T_w$ | Wall (surface) temperature |
| $T_\infty$ | Ambient temperature |

#### Convection (Mohammadkamal eq. 9)

$$q_{\text{conv}} = h \left( T_w - T_\infty \right)$$

where $h$ is the convective heat transfer coefficient [W/m²/K]. Set `hConv 0` to disable convection.

#### Evaporation (Proell 2023, eq. 12; Anisimov & Khokhlov 1995)

When the surface temperature exceeds the boiling point $T_v$, significant heat is lost through evaporative mass flux. To avoid numerical instability from the strong nonlinearity, the temperature used in the evaporation term is clipped: $[T] = \min(T_w, T_{\max})$ where $T_{\max} = T_v + 1000\,\text{K}$.

The evaporative mass flux is:

$$\dot{m} = 0.82 \, C_P \exp\left[ -C_T \left( \frac{1}{[T]} - \frac{1}{T_v} \right) \right] \sqrt{\frac{C_M}{[T]}}$$

and the evaporation heat flux is:

$$q_{\text{evap}} = \dot{m} \left( h_v + c_p \left([T] - T_{h,0}\right) \right), \quad [T] > T_v$$

| Symbol | Description | Typical value (Ti) |
|--------|-------------|-------------------|
| $C_P$ | Recoil pressure factor = $0.54 \, p_a$ | 54,700 Pa |
| $C_T$ | Temperature factor = $\bar{h}_v / R_{\text{gas}}$ | 56,600 K |
| $C_M$ | Factor $= M / (2\pi R_{\text{gas}})$ | 9.16×10⁻⁴ K·s²/m² |
| $h_v$ | Specific latent heat of evaporation | 9.83×10⁶ J/kg |
| $T_{h,0}$ | Enthalpy reference temperature | 663 K |
| $T_v$ | Boiling temperature (Ti) | 3560 K |
| $T_{\max}$ | Evaporation temperature cap | $T_v + 1000$ K |

The BC is implemented as a **pure Neumann** (`valueFraction = 0`) `mixedFvPatchField`:

```cpp
rGrad[faceI] = -qTot / max(kappa[faceI], SMALL);
valueFraction[faceI] = 0.0;
```

The material/heat-loss constants are read from
`constant/thermophysicalProperties`, not from the `T` field boundary entry.

#### Example `T` Boundary Entry

```
topWall
{
    type            pbfRadEvapTemperature;
    kappa           k;          // conductivity field name
    enableEvap      true;
    value           uniform 300;
}
```

#### Example `thermophysicalProperties` Entries

```
Tinf            Tinf [0 0 0 1 0 0 0] 300;
hConv           17;
epsilon         0.3;
Tv              3560;
Tmax            4560;
CP              54700;
CT              56600;
CM              9.16e-4;
hv              9.83e6;
Th0             663;
cpEvap          570;
```

### 6.2 Bottom Wall — Fixed Temperature (Dirichlet)

The bottom of the baseplate is held at the ambient temperature $T_\infty$ (Proell 2023, eq. 8):

$$T = T_\infty \quad \text{on } \Gamma_D$$

```
bottomWall { type fixedValue; value uniform 300; }
```

### 6.3 Side and Symmetry Walls — Thermally Insulating (Neumann)

All remaining boundaries are modeled as thermally insulating (zero heat flux), consistent with the assumption that the part is surrounded by low-conductivity powder (Proell 2023, eq. 9):

$$q \cdot n = 0 \quad \text{on } \Gamma_N$$

```
leftWall  { type zeroGradient; }
rightWall { type zeroGradient; }
front     { type zeroGradient; }
back      { type zeroGradient; }
```

---

## 7. Diffusion Number and Time Step Control

The solver monitors the **diffusion number** (analogous to the Fourier number) to assess the stability and accuracy of the time integration:

$$Di = \frac{\alpha \, \Delta t}{\Delta x^2}, \qquad \alpha = \frac{k}{\rho \, c_p^{\text{eff}}}$$

Although the implicit time scheme is unconditionally stable, a target diffusion number `maxDi` is enforced to control accuracy. The time step is adjusted adaptively in `controlDict`:

```
adjustTimeStep  yes;
maxDeltaT       1e-6;   // upper bound [s]
maxDi           0.5;    // target diffusion number
```

For scan-resolved LPBF simulations, the time step is also physically limited by the accuracy requirement that the laser should not travel more than its own radius in one time step (Proell 2023, eq. 20):

$$\Delta t \leq \frac{R}{v_{\text{scan}}}$$

This constraint holds regardless of the time integration scheme. For example, with $R = 50\,\mu\text{m}$ and $v_{\text{scan}} = 1\,\text{m/s}$, the maximum time step is $50\,\mu\text{s}$.

---

## 8. Solver Algorithm

Each time step executes the following sequence:

```
1. Store fL_old (previous liquid fraction)
2. updateThermo.H  → update ρ, cp, ks, km from temperature-dependent fits (Ti-6Al-4V)
3. updateLaser.H   → compute Q at new laser position and power
4. DiffusionNo.H   → compute Di; adjust Δt if adjustTimeStep=yes
5. Compute fL      → liquid fraction from T (linear ramp Ts→Tl)
6. Compute cpEff   → effective heat capacity with latent heat source term
7. Update rc       → consolidated fraction (max fL ever reached, irreversible)
8. Update rp, rm, rs → phase fractions from rc and fL
9. Update k        → phase-weighted conductivity
10. Update meltHistory → binary flag: 1 if cell ever reached Tl
11. Solve TEqn     → ρ cpEff ∂T/∂t = ∇·(k∇T) + Q
12. Write output fields (T, k, Q, rc, fL, rp, rm, rs, gradT, ...)
```

Step 11 uses OpenFOAM's `SIMPLE` non-orthogonal corrector loop to handle mesh non-orthogonality.

For multi-layer cases, the timestep loop is wrapped in an outer layer-activation loop driven by `constant/multiLayerProperties`. If `constant/dynamicMeshDict` is present, the base mesh can also be refined around the melt-history interface using `dynamicRefineFvMesh`.

---

## 9. Multi-Layer Activation and AMR

The multi-layer mode activates pre-meshed layers one at a time using `fvMeshSubset`. Layer 0 is seeded from the cells adjacent to `baseplatePatch`; subsequent layers are grown by face-neighbour expansion.

The key multi-layer inputs are read from `constant/multiLayerProperties`:

```text
baseplatePatch     baseplate;
nLayers            4;
cellsPerLayer      3;
layerDuration      200e-6;
layerInitialT      300;
exposedFacesPatch  top;
```

`cellsPerLayer` defaults to `1` if omitted. If `exposedFacesPatch` is not set, the exposed subset faces are written to the synthetic `oldInternalFaces` patch.

The solver keeps the persistent base-mesh fields under explicit names, interpolates them onto the active sub-mesh for each layer, and scatters the solved state back to the base mesh before writing. The multi-layer tutorial therefore writes the full domain at every output time, even though only the active subset is solved.

The AMR extension uses `dynamicRefineFvMesh` with a `refineIndicator` field built from smoothed `meltHistory`. A persistent `layerID` field tracks which layer each cell belongs to so layer membership survives mesh refinement.

The end-to-end example is `tutorials/multiLayer_basic`. Its `Allrun` script supports:

- `runMode=serial` for a single-process run
- `runMode=parallel` for the 4-rank setup used by the tutorial

---

## 10. Case Setup

### Directory Structure

```
case/
├── constant/
│   ├── LaserProperties          # laser model, rb, A, d, timeVsLaser* pointers
│   ├── thermophysicalProperties # material constants and poly coefficients
│   ├── timeVsLaserPosition      # table: (time   x y z)
│   └── timeVsLaserPower         # table: (time   Power[W])
├── system/
│   ├── blockMeshDict            # mesh definition
│   ├── controlDict              # time control, adjustTimeStep, maxDi
│   ├── fvSchemes                # discretization schemes
│   └── fvSolution               # linear solver settings
└── 0/ (or initial/)
    ├── T                        # initial temperature + boundary conditions
    ├── rc                       # initial consolidated fraction (0=powder, 1=solid)
    ├── rho                      # initial density
    ├── cp                       # initial specific heat
    └── Q                        # initial heat source (usually zeros)
```

### Laser Position and Power Tables

`constant/timeVsLaserPosition` defines the scan path as a time series:
```
(
    (0.0       0.0  0.0  0.001e-3)   // (time  x  y  z)
    (600e-6    6e-3 0.0  0.001e-3)
)
```

`constant/timeVsLaserPower` defines the power (supports pulsed-wave by alternating on/off):
```
(
    (0          195)    // laser on at 195 W
    (600e-6     195)
    (600.001e-6   0)    // laser off
)
```

### Initial Consolidated Fraction (`rc`)

- Set `rc = 1` for the solid substrate region.
- Set `rc = 0` for the powder layer.

This can be done with `setFields` using `setFieldsDict`.

### `LaserProperties` Dictionary

```
absorptivity    0.5;       // A: fraction of laser power absorbed
laserRadius     50e-6;     // R: beam radius [m]
absorptionDepth 35e-6;     // d: penetration depth [m]
laserModel      exponential;   // "cylindrical" or "exponential"

timeVsLaserPosition { file "$FOAM_CASE/constant/timeVsLaserPosition"; outOfBounds clamp; }
timeVsLaserPower    { file "$FOAM_CASE/constant/timeVsLaserPower";    outOfBounds clamp; }
```

### `thermophysicalProperties` Dictionary (Ti-6Al-4V example)

```
Ts    1877;       // solidus temperature [K]
Tl    1923;       // liquidus temperature [K]
L     2.85e5;     // latent heat of fusion [J/kg]
kp    0.2;        // powder conductivity [W/m/K]
ks    6.7;        // solid conductivity [W/m/K]  (initial; overwritten by updateThermo)
km    30.0;       // melt conductivity [W/m/K]   (initial; overwritten by updateThermo)
phi   0.35;       // powder porosity [-]
T1    1268;       // conductivity breakpoint [K]
// ... polynomial coefficients (see §5)
```

---

## 11. Building and Running

### Prerequisites

- OpenFOAM v2512 recommended by default for this repository
- The `pbfRadEvapTemperature` boundary condition library must be built first

### Build

```bash
# 1. Build the custom boundary condition library
cd applications/boundaryConditions/pbfRadEvapTemperature
wmake libso

# 2. Build the solver
cd ../../solvers/laserHeatFoam
wmake
```

### Run a Tutorial Case

```bash
cd tutorials/LPBF_titanium_case_195W_1mps
./Allrun
```

For the multi-layer tutorial:

```bash
cd tutorials/multiLayer_basic
./Allrun
```

or step by step:

```bash
blockMesh
setFields       # initialize rc for powder/substrate regions
laserHeatFoam   # run the solver
```

For parallel execution:

```bash
decomposePar
mpirun -np 4 laserHeatFoam -parallel
reconstructPar
```

---

## 12. Output Fields

The following fields are written every `writeInterval`:

| Field | Description |
|-------|-------------|
| `T` | Temperature [K] |
| `Q` | Volumetric laser heat source [W/m³] |
| `k` | Phase-weighted thermal conductivity [W/m/K] |
| `rc` | Consolidated fraction [0–1] (powder-to-solid history) |
| `fL` | Liquid fraction [0–1] |
| `rp` | Powder fraction [0–1] |
| `rm` | Melt fraction [0–1] |
| `rs` | Solid fraction [0–1] |
| `rho` | Effective density [kg/m³] |
| `cp` | Specific heat [J/kg/K] |
| `cpEff` | Effective specific heat including latent heat [J/kg/K] |
| `meltHistory` | Binary: 1 if cell ever exceeded `Tl` (full melt indicator) |
| `layerID` | Multi-layer membership index for each cell |
| `refineIndicator` | AMR driver field for `dynamicRefineFvMesh` |
| `gradT` | Temperature gradient vector [K/m] |

The `meltHistory`, `rc`, `layerID`, and `refineIndicator` fields are particularly useful for post-processing melt pool dimensions, track morphology, layer activation, and lack-of-fusion prediction.

---

## 13. References

```
[1] Mohammadkamal H., Caiazzo F. (2025).
    The role of laser operation mode on thermal and mechanical behavior
    in powder bed fusion: a numerical study.
    The International Journal of Advanced Manufacturing Technology, 140:3779–3796.
    https://doi.org/10.1007/s00170-025-16460-4

[2] Proell S.D., Munch P., Kronbichler M., Wall W.A., Meier C. (2023).
    A highly efficient computational approach for fast scan-resolved simulations
    of metal additive manufacturing processes on the scale of real parts.
    arXiv:2302.05164v3 [cs.CE].

[3] Proell S.D., Wall W.A., Meier C. (2020).
    On phase change and latent heat models in metal additive manufacturing
    process simulation.
    Advanced Modeling and Simulation in Engineering Sciences, 7:1–32.

[4] Anisimov S.I., Khokhlov V.A. (1995).
    Instabilities in laser-matter interaction. CRC Press.
    (evaporation model basis)
```

---

## License

This solver is based on OpenFOAM and is distributed under the **GNU General Public License v3**. See the OpenFOAM license for details.
