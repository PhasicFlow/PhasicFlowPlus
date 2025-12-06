# Unresolved Coupling System

## Table of Contents

1. [Overview](#1-overview)
2. [System Architecture](#2-system-architecture)
3. [Distribution Methods](#3-distribution-methods)
4. [Porosity Calculation](#4-porosity-calculation)
5. [Momentum Coupling](#5-momentum-coupling)
6. [Drag Force Models](#6-drag-force-models)
7. [Lift Force Models](#7-lift-force-models)
8. [Example Dictionary for Unresolved Coupling](#8-example-dictionary-for-unresolved-coupling)
9. [Nomenclature](#9-nomenclature)
10. [References](#10-references)

---

## 1. Overview

The unresolved coupling system module implements Eulerian-Lagrangian coupling for particle-fluid interactions in granular multiphase flows where the particle scale is not resolved in the fluid mesh. This approach is computationally efficient as it avoids the need for fine-scale fluid mesh around each particle, making it suitable for simulations with large numbers of particles.

## 2. System Architecture

### 2.1 Main Components

The unresolved coupling system consists of several key subsystems:

1. **Distribution Methods** - Maps particle properties to fluid cells using various smoothing techniques
2. **Porosity Calculation** - Computes fluid volume fraction in each cell
3. **Momentum Interactions** - Calculates drag, lift, and virtual mass forces
4. **Fluid Averaging (fluidVelocity)** - Evaluates fluid properties at particle locations
5. **Solid Averaging (solidVelocity)** - Averages particle properties over fluid cells

## 3. Distribution Methods

Distribution methods define how particle-scale data (velocity, force, etc.) are mapped to fluid cells and vice versa. Each method has different numerical characteristics and computational requirements.

### 3.1 Particle Centroid Method (PCM)

**Method Name:** `PCM`

The simplest distribution approach with no smoothing or distribution.

- Particle data is assigned directly to the cell containing the particle center
- No smoothing of data across adjacent cells
- Minimal computational cost
- Can produce localized fluctuations in cell-centered values
- Best suited for: dilute systems or when cell-to-particle size ratio is high (may be more than 6 or 7)

**Implementation:** When a particle center resides in a cell, all particle properties are assigned to that single cell without any distribution to neighboring cells.

### 3.2 Laplacian Diffusion Method

**Method Name:** `diffusion`

Uses Laplacian diffusion equation to smoothen particle data over multiple cells.

**Configuration Parameters:**
- `nSteps` - Number of diffusion steps to apply (integer)
- `standardDeviation` - Standard deviation parameter controlling diffusion extent (real)

**Mathematical Formulation:**

The diffusion equation smooths a scalar field $\phi$ as:

$$\frac{\partial \phi}{\partial \tau} = D \nabla^2 \phi \quad (1)$$

where:
- $D$ is the diffusion coefficient calculated as: $D = \sigma^2 / (4 \tau_{total})$ (see Eq. 2)
- $\sigma$ is the `standardDeviation` parameter
- $\tau_{total} = 1 \text{ s}$ is the total pseudo-time for integration
- $\tau$ is the pseudo-time variable for diffusion
- The equation is solved for `nSteps` iterations with time step $\Delta \tau = \tau_{total} / nSteps$ (see Eq. 3)

**Initial Condition:**

The initial field for the diffusion process is obtained using the **Particle Centroid Method (PCM)**, where particle properties are assigned directly to their containing cells without smoothing.

**Diffusion Coefficient:**

$$D = \frac{\sigma^2}{4 \tau_{total}} \quad (2)$$

**Time Step:**

$$\Delta \tau = \frac{\tau_{total}}{nSteps} \quad (3)$$

**Integration:**

The diffusion equation is integrated numerically up to pseudo-time $\tau = 1 \text{ s}$ using a finite volume scheme with `nSteps` iterations. This ensures consistent and reproducible smoothing independent of the grid size or computational parameters.

**Characteristics:**
- The recommended value for standard deviation is around 3 - 6 times particle diameter, depending on the cell-to-particle size ratio. 
- Provides smooth transitions of particle properties across cells
- Computationally more expensive than PCM due to multiple iterations
- Results in physically smoother field distributions
- Good balance between accuracy and computational cost
- Time-dependent diffusion ensures convergence to steady-state smooth distribution

**Sample Dictionary:**

```C++
diffusionInfo
{
    nSteps              5;
    standardDeviation   0.0075;
    
    // Optional: set to 1 to see solver output on screen
    log                 0;
}
```

**Parameters:**
- `nSteps` - Number of diffusion steps to apply (integer, required)
- `standardDeviation` - Standard deviation controlling diffusion extent (real, required)
- `log` - Optional parameter to enable solver output to screen (integer, default: 0; set to 1 to enable)

### 3.3 Gaussian Distribution Method

**Method Name:** `Gaussian`

Distributes particle data using a Gaussian (normal) distribution kernel centered at the particle position.

**Configuration Parameters:**
- `standardDeviation` - Standard deviation $\sigma$ of the Gaussian kernel (real)
- `maxLayers` - Maximum number of cell layers to consider (integer, default: 2)

**Mathematical Formulation:**

For a particle at position $\mathbf{x}_p$ with property value $q_p$, the weight assigned to cell $i$ is:

$$w_i = \frac{\exp\left(-\frac{|\mathbf{x}_p - \mathbf{x}_i|^2}{2\sigma^2}\right)}{\sum_j \exp\left(-\frac{|\mathbf{x}_p - \mathbf{x}_j|^2}{2\sigma^2}\right)} \quad (4)$$

The particle property contribution to cell $i$ is then:

$$q_i = w_i \cdot q_p \quad (5)$$

where $\sigma$ is controlled by the `standardDeviation` parameter.

**Averaged Property Calculation:**

The averaged property at the each cell center is calculated based on particles within a circle of radius $r = 3\sigma$ and center at the cell center using the Gaussian distribution:

$$\bar{q} = \sum_{p \in \text{circle}} q_p w_p \quad (5a)$$

where $w_p$ is defined in Eq. (4) and normalized such that $\sum_p w_p = 1$ within the search circle. The search radius $r = 3\sigma$ ensures that only particles with significant Gaussian weight contribution are included in the averaging, reducing computational cost while maintaining accuracy.

**Characteristics:**
- Smooth Gaussian-shaped distribution kernel
- Weights decrease exponentially with distance
- Requires neighbor cell list construction
- Good for capturing smooth spatial gradients
- The `maxLayers` parameter defines the search radius for neighboring cells

**Sample Dictionary:**

```C++
GaussianInfo
{
    standardDeviation   0.0075;
    maxLayers           2;
}
```

**Parameters:**
- `standardDeviation` - Standard deviation $\sigma$ of the Gaussian kernel (real, required)
- `maxLayers` - Maximum number of cell layers to consider for neighbor search (integer, default: 2)

### 3.4 Gaussian Integral Method

**Method Name:** `GaussianIntegral`

Uses Gaussian integral for distributing/smoothing particle data, providing more accurate volume-weighted distribution than simple point evaluation.

**Configuration Parameters:**
- `maxLayers` - Maximum number of cell layers to consider (integer, default: 1)

**Characteristics:**
- Integrates Gaussian distribution over cell volumes
- More numerically accurate than point-based Gaussian distribution
- Higher computational cost due to integration
- Better conservation properties
- Recommended for high-accuracy requirements

**Reference:** Kianimoqadam, A., Lapp, J. (2026). "Gaussian integral method for void fraction." Particuology, 108, 125-142. [https://doi.org/10.1016/j.partic.2025.10.014](https://doi.org/10.1016/j.partic.2025.10.014)

**Sample Dictionary:**

```C++
GaussianIntegralInfo
{
    // optional, default is 1
    maxLayers           2;
}
```

### 3.5 Adaptive Gaussian Method

**Method Name:** `adaptiveGaussian`

Uses Gaussian distribution with adaptive standard deviation that varies based on the particle-to-cell size ratio.

**Mathematical Formulation:**

The characteristic cell size is calculated as:

$$d_{cell} = V_{cell}^{1/3} \quad (6)$$

The standard deviation is adaptively calculated as:

$$\sigma = d_{cell} \cdot f_s \cdot a \cdot \left(\frac{d_{cell}}{d_p}\right)^{e} \quad (7)$$

where:
- $d_{cell}$ is the characteristic cell size (Eq. 6)
- $f_s$ is a smoothing factor (default: 1.0, optional in configuration)
- $a = 0.6142275$ (empirical constant)
- $e = -0.6195039$ (empirical exponent)
- $d_p$ is the particle diameter

After calculating the adaptive standard deviation, the method applies Eq. (4) and Eq. (5a) with the adjusted $\sigma$ value to compute the distribution weights and averaged properties.

**Characteristics:**
- The method is accurate for cell-to-particle size ratios down to 1. It can be used in almost every unresolved simulation.
- In the regions where mesh non-orthogonality increases, the error of this method increaases. 
- Automatically adjusts smoothing based on size distribution
- Useful for polydisperse systems with wide particle size ranges
- Provides consistent behavior across different Cell sizes
- Efficient in computations compared to `diffusion` and `subDivision29`. 


**Sample Dictionary:**

```C++
adaptiveGaussianInfo
{    
    maxLayers           1;          // optional default: 1

    smoothingFactor     1.0;        // optional, default: 1.0
}
```

### 3.6 Sub-Division Methods

#### 3.6.1 29-Sub-volume Method

**Method Name:** `subDivision29`

Divides a sphere into 29 equal parts and locates each part in cells.

**Characteristics:**

- Divides particle sphere into 29 equal volumetric segments
- Locates each segment in the appropriate cell based on intersection
- Accounts for partial particle volumes in cells
- High accuracy for porosity calculations
- Higher computational cost due to geometric calculations
- Best suited for: Accurate porosity calculations in dense suspensions with cell-to-particle size ration higher than 3.

#### 3.6.2 9-Sub-volumne Method

**Method Name:** `subDivision9`

Divides a sphere into 9 equal parts.

**Characteristics:**

- Simpler version with 9 segments instead of 29.
- Trade-off between accuracy and computational efficiency
- Suitable for moderate accuracy requirements
- Lower computational cost than 29-subdivision
- Use it for simulations with cell-to-particle size ratio higher than 4.

## 4. Porosity Calculation

Porosity (fluid volume fraction) is a critical parameter in unresolved coupling. It determines the local volume of fluid available in each cell.

### 4.1 Available Porosity Methods

#### 4.1.1 Porosity with Distribution (Default)

**Method Name:** `distribution`

Uses the selected distribution method to map particle volume to cells.

**Configuration:**
```
porosity
{
    method      distribution;
    alphaMin    0.2;  // Minimum allowed porosity
}
```

**Calculation:**
For each particle with diameter $d_p$:

**Particle volume:**

$$V_p = \frac{\pi}{6} d_p^3 \quad (8a)$$

**Fluid volume fraction in cell:**

$$\alpha = 1 - \frac{V_{solid,cell}}{V_{cell}} \quad (8b)$$

where $V_{solid,cell}$ is the total solid volume in the cell (distributed from particles using the selected distribution method), and $V_{cell}$ is the cell volume.

**Calculation of $V_{solid,cell}$:**

The solid volume in each cell is calculated by summing the weighted contributions from all particles:

$$V_{solid,cell} = \sum_{p} w_{p,cell} \cdot V_p \quad (8c)$$

where:
- $w_{p,cell}$ is the distribution weight for particle $p$ in cell (determined by the selected distribution method, e.g., PCM, Diffusion, Gaussian, etc.)
- $V_p$ is the particle volume (Eq. 8a)
- The sum is over all particles that contribute to the cell

The distribution weights ensure that:
- **PCM**: All particle volume goes to the cell containing the particle center ($w_{p,cell} = 1$ for target cell, 0 elsewhere)
- **Diffusion**: Volume is spread using Laplacian diffusion weights (Eq. 2-3)
- **Gaussian**: Volume is distributed using Gaussian kernel weights (Eq. 4-5a)
- **GaussianIntegral**: Volume is distributed using volume-weighted Gaussian integration
- **AdaptiveGaussian**: Volume is distributed using adaptively adjusted Gaussian weights (Eq. 7)

The minimum porosity `alphaMin` is enforced to prevent numerical issues:

$$\alpha_{final} = \max(\alpha_{calculated}, \alpha_{min}) \quad (9)$$

#### 4.1.2 Sub-Division 29

**Method Name:** `subDivision29`

Uses the 29-sub-division method for porosity calculation with high accuracy.

#### 4.1.3 Sub-Division 9

**Method Name:** `subDivision9`

Uses the 9-sub-division method for porosity calculation with good balance of accuracy and speed.

## 5. Momentum Coupling

Momentum coupling transfers forces between particles and fluid and calculates interaction forces.

**General Momentum Source Term:**

The momentum source term applied to the fluid phase is expressed in the form:

$$S = S_p \mathbf{U} + S_u \quad (10)$$

where:
- $S_p$ is the implicit coefficient (related to drag and other velocity-dependent forces)
- $\mathbf{U}$ is the fluid velocity
- $S_u$ is the explicit source term

The coefficients $S_p$ and $S_u$ are computed when drag force is calculated:

$$F_D = \frac{V_p \beta}{1-\alpha}(\mathbf{U}-v_p)$$

The momentum source in each cell is calculated by summing the weighted contributions from all particles:

$$S_p = \frac{1}{V_c} \sum_{p} w_{p,c} \frac{V_p \beta}{1-\alpha} \quad (11)$$

$$S_u = - \frac{1}{V_c} \sum_{p} w_{p,c} \frac{V_p \beta}{1-\alpha} v_p \quad (12)$$

In the calculations of other forces, they are added to the explicit part of the momentum source term:

$$S_u = \frac{1}{V_c} \sum_{p} w_{p,c} \left( - \frac{V_p \beta}{1-\alpha} v_p + \mathbf{F}_{virtual,p} + \mathbf{F}_{lift,p} \right) \quad (13)$$

where:
- $w_{p,c}$ is the distribution weight for particle $p$ in cell $c$ (from the selected distribution method, Eq. 4)
- $V_p$ is the particle volume
- $V_c$ is the cell volume
- $\beta$ is the drag friction (from selected drag model in Section 6)
- $\alpha$ is the local porosity
- $v_p$ is the particle velocity
- $\mathbf{F}_{virtual,p}$ is the virtual mass force on particle $p$ (inertia of displaced fluid)
- $\mathbf{F}_{lift,p}$ is the lift force on particle $p$ (from selected lift model in Section 7)
- The sums are over all particles that contribute to the cell

### 5.1 Momentum Exchange Mapping

Controls how momentum exchange terms are distributed to fluid cells (based on Eqs. (11-13)).

**Options:**

- **`cell`** - No smoothing; momentum exchange assigned to cell containing particle center (Uses PCM).
- **`distribution`** - Uses the selected distribution method for smoothing momentum terms across cells.

### 5.2 Fluid Velocity Evaluation

Determines how fluid velocity at particle center is evaluated.

**Options:**

- **`cell`** - Uses cell value (the cell containing particle center)
- **`interpolate`** - Interpolates from surrounding cell values to particle center
- **`distribution`** - Uses distribution method to obtain volume-averaged velocity at particle center

**Interpolate Method:**

The fluid velocity at particle position is computed using inverse distance weighting from the particle cell and its neighboring cells:

$$\overline{\mathbf{U}}_p = \frac{\sum_{c} \frac{1}{r_{pc}} \mathbf{U}_c}{\sum_{c} \frac{1}{r_{pc}}} \quad (14a)$$

where:
- $r_{pc}$ is the distance between particle $p$ and cell center $c$
- $\mathbf{U}_c$ is the fluid velocity at cell center $c$
- The sum includes the particle's containing cell and all neighboring cells
- $r_{pc}$ is bounded below by a small value (SMALL) to avoid division by zero

**Distribution Method:**

The fluid velocity at the particle center is obtained by applying the inverse distribution operation (weighted averaging) using the distribution weights:

$$\overline{\mathbf{U}}_p = \sum_{c} w_{p,c} \mathbf{U}_c \quad (14b)$$

where:
- $w_{p,c}$ is the distribution weight for particle $p$ in cell $c$ (from the selected distribution method, Eq. 4)
- $\mathbf{U}_c$ is the fluid velocity at cell center $c$
- The sum is over all cells that contribute to the particle (determined by the distribution method)
- This approach is consistent with how particle properties are distributed to cells, providing a symmetric averaging scheme

### 5.3 Solid Velocity Evaluation

Determines how particle velocity is evaluated in coupling calculations.

**Options:**

- **`particle`** - Uses exact particle velocity
- **`distribution`** - Uses distribution method to obtain volume-averaged velocity in the cell; this averaged value is used for all particles in the cell

**Particle Method:**

The particle velocity is used directly without averaging:

$$\overline{\mathbf{v}}_c = \mathbf{v}_p \quad (15a)$$

where:
- $\mathbf{v}_p$ is the exact velocity of particle $p$
- This is applied when calculating momentum coupling for a particle in cell $c$

**Distribution Method:**

The velocity is first distributed and accumulated in cells, then normalized by the solid volume fraction in the cell:

$$\overline{\mathbf{v}}_c = \frac{\sum_p V_p w_{p,c} \mathbf{v}_p}{(1-\alpha_c) V_c} \quad (15b)$$

Then this cell-averaged value is used for all particles in that cell:

$$\overline{\mathbf{v}}_{p} = \overline{\mathbf{v}}_c \quad \text{for all particles } p \text{ in cell } c \quad (15c)$$

where:
- $V_p$ is the particle volume
- $w_{p,c}$ is the distribution weight for particle $p$ in cell $c$ (from Eq. 4)
- $\mathbf{v}_p$ is the velocity of particle $p$
- $\alpha_c$ is the porosity in cell $c$
- $V_c$ is the cell volume
- The first equation accumulates volume-weighted particle velocities, normalized by solid volume in the cell
- The second equation shows that this averaged velocity is then applied to all particles in the cell for momentum coupling calculations

## 6. Drag Force Models

Drag force is the primary interaction force in unresolved coupling, representing resistance to relative motion between particles and fluid.

### 6.1 Drag Force Parameters

Two key parameters are required for all drag force calculations:

**Relative Velocity:**

The relative velocity between fluid and particle is obtained from the averaged velocities determined in Section 5:

$$\mathbf{U}_{rel} = \overline{\mathbf{U}}_f - \overline{\mathbf{v}}_p \quad (16a)$$

where:
- $\overline{\mathbf{U}}_f$ is the fluid velocity evaluated at the particle position (using interpolate or distribution method, Eq. 14a or 14b)
- $\overline{\mathbf{v}}_p$ is the particle velocity (either particle method giving $\mathbf{v}_p$ directly, or distribution method giving $\overline{\mathbf{v}}_c$, Eq. 15a or 15c)

**Reynolds Number:**

The particle Reynolds number is calculated from the relative velocity magnitude, particle diameter, and fluid kinematic viscosity:

$$Re = \alpha \frac{\rho_f |\mathbf{U}_{rel}| d_p}{\mu_f} \quad (16b)$$

where:
- $\alpha$ is the local porosity (fluid volume fraction) in the cell containing the particle
- $\rho_f$ is the fluid density
- $|\mathbf{U}_{rel}|$ is the magnitude of the relative velocity
- $d_p$ is the particle diameter (from Eq. 8a)
- $\mu_f$ is the fluid dynamic viscosity

The porosity factor $\alpha$ accounts for the concentration effects on drag in the unresolved Eulerian-Lagrangian coupling framework. A minimum threshold value `residualRe` is applied in implementations to prevent numerical issues at very low Reynolds numbers: $Re \geq Re_{min}$.

**Relationship between Dimensionless Drag and Drag Coefficient:**

The dimensionless drag function $\hat{f}^d$ computed by each drag model is related to the drag coefficient $\beta$ used in momentum coupling through:

$$18\mu_f \alpha (1-\alpha) \hat{f}^d(\alpha, Re) = \beta d_p^2 \quad (16c)$$

Solving for the drag coefficient:

$$\beta = \frac{18\mu_f \alpha (1-\alpha)}{d_p^2} \hat{f}^d(\alpha, Re) \quad (16d)$$

where:
- $\hat{f}^d(\alpha, Re)$ is the dimensionless drag function from the drag model correlation
- $\beta$ is the drag coefficient appearing in momentum source terms (Eq. 11-13)
- $\mu_f$ is the fluid dynamic viscosity
- $\alpha$ is the fluid volume fraction (porosity)
- $(1-\alpha)$ is the solid volume fraction
- $d_p$ is the particle diameter

The drag force on each particle is then:


$$F_D = \frac{V_p \beta}{1-\alpha}(\mathbf{U}-v_p) \quad (16e)$$

where $V_p$ is the particle volume.

### 6.2 Drag Model Closures

All drag models accept:
- `model` - Selection of drag closure (e.g., `DiFelice`, `ErgunWenYu`, `Beetstra`, `Rong`, `Cello`)
- `residualRe` - Minimum Reynolds number threshold to prevent numerical issues at very low $Re$

### 6.3 DiFelice Correlation

**Model Name:** `DiFelice`

Di Felice (1994) provides a correlation for the dimensionless drag function $\hat{f}^d$ accounting for void fraction effects.

**Formulation:**

$$\hat{f}^d = \frac{C_d}{24} Re \cdot \alpha^{-\xi} \quad (17)$$

where:
- $\hat{f}^d$ is the **dimensionless drag function** computed from this correlation
- $C_d = \left(0.63 + \frac{4.8}{\sqrt{Re}}\right)^2$ is the shape drag coefficient
- $Re$ is the particle Reynolds number (Eq. 16b)
- $\alpha$ is the porosity (fluid volume fraction)
- $\xi = 3.7 - 0.65 \exp\left(-0.5(1.5 - \log_{10}Re)^2\right)$ is the voidage exponent

The drag coefficient $\beta$ is obtained from Eq. 16d using this dimensionless drag value.

**Characteristics:**
- Widely used and well-validated
- Good for low to intermediate Reynolds numbers (should be of primary choice for general applications)
- Accounts for local porosity effects through $\alpha^{-\xi}$

**Reference:** Di Felice, R. (1994) "The voidage function for fluid–particle interaction systems." International Journal of Multiphase Flow, 20, 153–159.

### 6.4 Ergun-Wen-Yu Correlation

**Model Name:** `ErgunWenYu`

Combines Ergun equation for dense flows with Wen-Yu correlation for dilute flows to provide the dimensionless drag function $\hat{f}^d$.

**Formulation:**

For dilute flows ($\alpha \geq 0.8$):

$$\hat{f}^d = \frac{C_d}{24} Re \alpha^{-3.65} \quad (18)$$

where:
$$C_d = \begin{cases}
24(1 + 0.15 Re^{0.687})/Re & \text{if } Re \leq 1000 \\
0.44 & \text{if } Re > 1000
\end{cases}$$

For dense flows ($\alpha < 0.8$):

$$\hat{f}^d = \frac{150(1-\alpha)}{18 \alpha^2} + \frac{1.75 Re}{18 \alpha^2} \quad (19)$$

where $\hat{f}^d$ is the **dimensionless drag function** for the respective flow regime.


**Characteristics:**
- Handles both dense and dilute regimes
- Transitions smoothly at $\alpha = 0.8$
- Particularly suitable for packed bed simulations as well as fluidized beds

**References:**
- Gidaspow, D. (1994) "Multiphase Flow and Fluidization"
- Ergun, S. (1952) "Fluid flow through packed columns"
- Wen, C.Y. and Yu, Y.H. (1966) "Mechanics of fluidization"

### 6.5 Beetstra Correlation

**Model Name:** `Beetstra`

Drag correlation providing the dimensionless drag function $\hat{f}^d$ for intermediate Reynolds number flows past mono-disperse arrays of spheres.

**Formulation:**

$$\hat{f}^d = \frac{10(1-\alpha)}{\alpha^2} + \alpha^2(1+1.5\sqrt{1-\alpha}) + \frac{0.413}{24\alpha^2} \cdot \frac{\left(\frac{1}{\alpha} + 3\alpha(1-\alpha) + 8.4 Re^{-0.343}\right)}{1 + 10^{3(1-\alpha)} Re^{-0.5 + 2(1-\alpha)}} \cdot Re \quad (20)$$

where $\hat{f}^d$ is the **dimensionless drag function** from this complex correlation.

**Characteristics:**
- High accuracy across range of porosities and Reynolds numbers
- More complex formula with better accuracy than simpler models
- Computationally more expensive

**Reference:** Beetstra, R., van der Hoef, M. A., & Kuipers, J. A. M. (2007). "Drag force of intermediate Reynolds number flow past mono- and bidisperse arrays of spheres." AIChE Journal, 53(3), 489-501.

### 6.6 Rong Correlation

**Model Name:** `Rong`

Drag correlation based on lattice-Boltzmann simulations providing the dimensionless drag function $\hat{f}^d$ for uniform sphere packing.

**Formulation:**

$$\hat{f}^d = \frac{C_d}{24} Re \cdot \alpha^{-\xi} \quad (21)$$

where:
- $\hat{f}^d$ is the **dimensionless drag function** from lattice-Boltzmann data
- $C_d = (0.63 + 4.8/\sqrt{Re})^2$ is the shape drag coefficient
- $\xi = 2.65(1+\alpha) - (5.3-3.5\alpha)\alpha^2\exp\left(-0.5(1.5-\log_{10}(Re))^2\right)$ is the porosity-dependent exponent

The drag coefficient $\beta$ is obtained from Eq. 16d using this dimensionless drag value.

**Characteristics:**
- Based on detailed LB simulations
- Porosity-dependent exponent accounts for packing effects

**Reference:** Rong, L.W., Dong, K.J., Yu, A.B. (2013) "Lattice-Boltzmann simulation of fluid flow through packed beds of uniform spheres." Chemical Engineering Science, 99, 44-58.

## 7. Lift Force Models

Lift forces arise from relative motion and rotation effects and can be significant in certain flow regimes.

**General Form:**

$$\mathbf{F}_L = \frac{\pi d^2}{8} \rho_f |\mathbf{U}_{rel}|^2 (C_{L,shear} \mathbf{e}_{shear} + C_{L,spin} \mathbf{e}_{spin}) \quad (22)$$

### 7.1 Lift Force Parameters

The general lift force form involves several key parameters that must be computed from fluid and particle properties. The definitions below are derived directly from the source code implementations.

**Relative Velocity:**

The relative velocity between the particle and the fluid is computed as:

$$\mathbf{U}_{rel} = \mathbf{v}_p - \mathbf{U}_f \quad (22a)$$

where:
- $\mathbf{v}_p$ is the particle velocity
- $\mathbf{U}_f$ is the fluid velocity at the particle location

Note: This differs from Eq. 16a (used in drag calculations) by a sign; here the direction is particle velocity minus fluid velocity, which determines the reference frame for lift calculations.

**Fluid Vorticity:**

The fluid vorticity is computed as the curl of the velocity field:

$$\boldsymbol{\omega}_f = \nabla \times \mathbf{U}_f \quad (22b)$$

This represents the local rotation (vorticity) of the fluid, which is a key parameter for understanding shear-induced lift.

**Shear-Induced Lift Direction:**

The unit vector in the direction of shear-induced lift is determined by the cross product of the fluid vorticity and the relative velocity vector:

$$\mathbf{e}_{shear} = \frac{\boldsymbol{\omega}_f \times \mathbf{U}_{rel}}{|\boldsymbol{\omega}_f \times \mathbf{U}_{rel}|} \quad (22c)$$

This direction represents the Magnus effect due to the ambient fluid rotation. It is perpendicular to both the vorticity vector and the relative velocity vector.

**Spin-Induced Lift Direction:**

The unit vector in the direction of spin-induced (Magnus) lift is determined by the cross product of the particle angular velocity and the relative velocity vector:

$$\mathbf{e}_{spin} = \frac{\boldsymbol{\omega}_p \times \mathbf{U}_{rel}}{|\boldsymbol{\omega}_p \times \mathbf{U}_{rel}|} \quad (22d)$$

where $\boldsymbol{\omega}_p$ is the particle angular velocity (rotational velocity). This direction represents the Magnus effect due to particle rotation in the relative flow field.

**Reynolds Numbers and Dimensionless Parameters:**

The following Reynolds numbers and dimensionless parameters are used across all lift models:

- $Re_p = \frac{|\mathbf{U}_{rel}| d_p}{\nu}$ - Particle Reynolds number
- $Re_\omega = \frac{|\boldsymbol{\omega}_p| d_p^2}{\nu}$ - Particle rotation Reynolds number  
- $Re_\Omega = \frac{|\nabla \times \mathbf{U}_f| d_p^2}{\nu}$ - Fluid vorticity Reynolds number
- $Sr = \frac{Re_\Omega}{Re_p}$ - Shear rate parameter
- $Rr = \frac{Re_\omega}{Re_p}$ - Rotation rate parameter
- $\epsilon = \frac{\sqrt{Re_\omega}}{Re_p}$ - Dimensionless rotation parameter

where $d_p$ is the particle diameter and $\nu$ is the kinematic viscosity of the fluid.

### 7.2 Saffmann Lift Model

**Model Name:** `Saffmann`

Low Reynolds number lift model based on Saffman (1965).

**Applicability:**
- Particle Reynolds number: $Re_p < 1$
- Low shear rate conditions

**Formulation:**

Shear lift coefficient:

$$C_{L,shear} = \frac{18}{\pi^2}\sqrt{\frac{Sr}{Re_p}} J - \frac{11}{8} Sr \exp(-0.5 Re_p) \quad (23)$$

where $J = 2.255$.

Spin lift coefficient:

$$C_{L,spin} = Rr \quad (24)$$

**Reference:** Saffman, P.G.T., 1965. "The lift on a small sphere in a slow shear flow." J. Fluid Mech. 22, 385–400.

### 7.3 Loth2008 Lift Model

**Model Name:** `Loth2008`

Wide-range lift model accounting for both shear-induced and spin-induced effects with capability for negative lift prediction.

**Spin-induced Lift Coefficient:**

$$C_{L\Omega} = Rr \left( 1 - \{0.675 + 0.15(1 + \tanh[0.28(Rr-2)])\} \tanh(0.18 \sqrt{Re_p}) \right) \quad (25)$$

**Shear-induced Lift Coefficient:**

For $Re_p \leq 50$:

$$C_{L\omega} = \frac{18}{\pi^2} \sqrt{\frac{Sr}{Re_p}} J(\epsilon) \quad (26)$$

where:
$$J(\epsilon) = \begin{cases}
-0.04\epsilon + 2.05\epsilon^2 - 32.2\epsilon^3 + 106.8\epsilon^4 & \text{if } \epsilon \leq 0.23 \\
\frac{2.225}{(1 + 0.02304/\epsilon^2)^{12.77}} & \text{if } \epsilon > 0.23
\end{cases} \quad (27)$$

For $Re_p > 50$:

$$C_{L\omega} = -Sr^{1/3}\left(0.0525 + 0.0575 \tanh(5 \log_{10}(Re_p/120))\right) \quad (28)$$

**Characteristics:**
- Covers wide range of Reynolds numbers
- Accurately predicts negative lift in certain regimes
- More complex than Saffmann model but more generally applicable

**Reference:** Loth, E., 2008. "Lift of a spherical particle subject to vorticity and/or spin." AIAA Journal, 46(4), 801–809.

### 7.4 Shi2019 Lift Model

**Model Name:** `Shi2019`

Advanced lift force model based on Shi and Rzehak (2019) for wide-range Reynolds number applicability.

**Spin Lift Coefficient:**

$$C_{L,spin} = Rr \left( 1 - 0.62 \tanh(0.3 \sqrt{Re_p}) - 0.24 \frac{\tanh(0.01 Re_p)}{\tanh(0.8 \sqrt{Rr})} \arctan[0.47(Rr-1)] \right) \quad (29)$$

**Fluid Shear Lift Coefficient (for $Re_p \leq 50$):**

With $\epsilon = \sqrt{Re_\omega} / Re_p$:

If $\epsilon \leq 0.23$:
$$J = -0.04\epsilon + 2.05\epsilon^2 - 32.2\epsilon^3 + 106.8\epsilon^4 \quad (30)$$

If $\epsilon > 0.23$:
$$J = \frac{2.225}{(1 + 0.02304/\epsilon^2)^{12.77}} \quad (31)$$

$$C_{L,shear} = \frac{18}{\pi^2} \sqrt{\frac{Sr}{Re_p}} J(\epsilon) - \frac{11}{8} Sr \exp(-0.5 Re_p) \quad (32)$$

**Fluid Shear Lift Coefficient (for $Re_p > 50$):**

$$C_{L,shear} = -0.064 \exp(0.525 Sr) \left( 0.49 + 0.51 \tanh \left[ 5\log_{10}\left(\frac{Re_p Sr^{0.08}}{120}\right) \right] \right) \quad (33)$$

**Characteristics:**
- Latest lift model with advanced flow regime handling
- Improved shear lift prediction at high Reynolds numbers
- Accounts for complex interactions between shear and spin
- Recommended for: High-Reynolds number systems with significant shear

**Reference:** Pengyu Shi and Roland Rzehak, 2019. "Lift forces on solid spherical particles in unbounded flows." Chemical Engineering Science, 208, 115145.

### 7.5 Surface Rotation Torque Models

Models for torque on particles due to surface rotation effects.

**Available Options:**
- **`none`** - Not included (default)
- **`lowReynolds`** - Based on Happel and Brenner model for low Reynolds numbers
- **`Loth2008`** - Based on Loth (2008) model
- **`Shi2019`** - Based on Shi model



## 8. Example Dictionary for Unresolved Coupling

A complete configuration example (sampleDictionary):

```C++
unresolved
{
    // Distribution method for mapping particle data (porosity, velocity,
    // force, etc.) over fluid cells
    // 
    // Available methods: 
    //     - PCM: Particle Centroid Method (no smoothing)
    //     - diffusion: Laplacian diffusion for smoothing
    //     - Gaussian: Gaussian distribution with specified std dev
    //     - GaussianIntegral: Gaussian integral for distributing data
    //     - AdaptiveGaussian: Gaussian with adaptive std deviation
    //     - subDivision29: Divided-volume method (29-sub-volume version)
    //     - subDivision9: Divided-volume method (9-sub-volume version)
    distributionMethod      diffusion;
    
    // Distribution method required settings 
    diffusionInfo
    {
        nSteps              5;
        standardDeviation   0.0075;
    }

    // Required settings for calculating porosity method 
    porosity
    {
        // method is optional
        //    - default value is distribution
        //    - Other options: subDivision29, subDivision9 
        method      distribution; 

        // alphaMin is minimum alpha allowed in porosity calculations
        alphaMin    0.2;
    }

    // Settings for momentum coupling  
    momentumInteraction
    {
        // How to perform smoothing on momentum exchange terms
        // Available options are: 
        //    - cell: No smoothing, assigned to containing cell
        //    - distribution: Uses distributionMethod for smoothing
        momentumExchange distribution; 

        // How to evaluate fluid velocity at particle location
        // Available options are: 
        //    - cell: Uses cell value
        //    - interpolate: Interpolates from neighboring cells
        //    - distribution: Uses distributionMethod for averaging
        fluidVelocity    distribution;

        // How to evaluate particle velocity in calculations
        // Available options are:
        //    - particle: Uses exact particle velocity
        //    - distribution: Uses distributionMethod to obtain average
        //      velocity in the cell for all particles
        solidVelocity    distribution;

        drag
        {
            // Drag force closure
            // Available options for spherical particles:
            //   No size distribution treatment:
            //     - DiFelice, ErgunWenYu, Rong, Cello, Beetstra
            //   For systems with size distribution:
            //     - CelloPolydisperse, BeetstraPolydisperse, adhoc
            model       DiFelice; 

            // Residual Reynolds number 
            residualRe  1.00e-6;
        }

        lift
        {   
            // Lift force model options:
            //   - none: Not included (default)
            //   - Saffmann: Saffman (1965), low Re
            //   - Loth2008: Loth (2008), intermediate Re
            //   - Shi2019: Shi & Rzehak (2019), wide Re range
            // References:
            //   Saffman, P.G.T., 1965. J. Fluid Mech. 22, 385-400
            //   Loth, E., 2008. AIAA J. 46, 801-809
            //   Shi, P. & Rzehak, R., 2019. Chem. Eng. Sci. 208,
            //     115145
            model                   Loth2008;

            // Surface rotation torque model options:
            //    - none: Not included (default)
            //    - lowReynolds: Happel and Brenner model
            //    - Loth2008: Based on Loth (2008) model
            //    - Shi2019: Based on Shi model
            surfaceRotationTorque   lowReynolds; 

            residualRe              1.0e-6;
            
        }

        virtualMass
        {
            // This part has not been implemented yet
        }
    }

    heatTransfer
    {

    }

    massTransfer
    {}
    
}
```

## 9. Nomenclature

### Greek Symbols

| Symbol | Definition | Units |
|--------|-----------|-------|
| $\alpha$ | Fluid volume fraction (porosity) | - |
| $\alpha_{min}$ | Minimum allowed porosity | - |
| $\beta$ | Interphase momentum transfer/fluid friction coefficient | kg/(m³·s) |
| $\Delta \tau$ | Time step in pseudo-time integration | s |
| $\epsilon$ | Dimensionless rotation parameter = $\sqrt{Re_\omega}/Re_p$ | - |
| $\xi$ | Voidage exponent in DiFelice correlation | - |
| $\mu_f$ | Dynamic viscosity of fluid | Pa·s |
| $\nu$ | Kinematic viscosity of fluid = $\mu_f/\rho_f$ | m²/s |
| $\rho_f$ | Density of fluid | kg/m³ |
| $\rho_p$ | Density of particle | kg/m³ |
| $\sigma$ | Standard deviation of Gaussian distribution kernel | m |
| $\tau$ | Pseudo-time variable in diffusion | s |
| $\tau_{total}$ | Total pseudo-time for diffusion integration | s |
| $\boldsymbol{\omega}_f$ | Fluid vorticity vector = $\nabla \times \mathbf{U}_f$ | rad/s |
| $\boldsymbol{\omega}_p$ | Particle angular velocity (rotational velocity) | rad/s |

### Roman Symbols (Scalars)

| Symbol | Definition | Units |
|--------|-----------|-------|
| $a$ | Empirical constant in adaptive Gaussian method | - |
| $C_d$ | Drag coefficient (shape-dependent) | - |
| $C_{L,shear}$ | Shear-induced lift coefficient | - |
| $C_{L,spin}$ | Spin-induced (Magnus) lift coefficient | - |
| $C_{L\Omega}$ | Spin-induced lift coefficient (Loth2008 notation) | - |
| $C_{L\omega}$ | Shear-induced lift coefficient (Loth2008 notation) | - |
| $d_c$ | Characteristic cell size = $V_{cell}^{1/3}$ | m |
| $d_p$ | Particle diameter | m |
| $D$ | Diffusion coefficient | m²/s |
| $e$ | Empirical exponent in adaptive Gaussian method | - |
| $f_s$ | Smoothing factor in adaptive Gaussian method | - |
| $\hat{f}^d$ | Dimensionless drag function | - |
| $F_D$ | Drag force magnitude on particle | N |
| $J$ | Function in Saffmann lift correlation | - |
| $J(\epsilon)$ | Function of dimensionless rotation parameter | - |
| $nSteps$ | Number of diffusion or integration steps | - |
| $r_{pc}$ | Distance between particle $p$ and cell center $c$ | m |
| $Re$ | Particle Reynolds number | - |
| $Re_p$ | Particle Reynolds number in lift calculations | - |
| $Re_s$ | Shear Reynolds number (fluid vorticity) | - |
| $Re_w$ | Rotation Reynolds number (particle spin) | - |
| $Re_\omega$ | Rotation Reynolds number (particle angular velocity) | - |
| $Re_\Omega$ | Fluid vorticity Reynolds number | - |
| $Rr$ | Rotation rate parameter = $Re_\omega/Re_p$ | - |
| $S$ | Momentum source term | N/m³ |
| $S_p$ | Implicit coefficient in momentum source | kg/m³ |
| $S_u$ | Explicit part of momentum source | N/m³ |
| $Sr$ | Shear rate parameter = $Re_\Omega/Re_p$ | - |
| $V_c$ | Cell volume | m³ |
| $V_{cell}$ | Cell volume | m³ |
| $V_p$ | Particle volume = $\pi d_p^3/6$ | m³ |
| $V_{solid,cell}$ | Total solid volume in cell | m³ |

### Roman Symbols (Vectors)

| Symbol | Definition | Units |
|--------|-----------|-------|
| $\mathbf{e}_{shear}$ | Unit vector in direction of shear-induced lift | - |
| $\mathbf{e}_{spin}$ | Unit vector in direction of spin-induced lift | - |
| $\mathbf{F}_{lift,p}$ | Lift force vector on particle $p$ | N |
| $\mathbf{F}_{virtual,p}$ | Virtual mass force on particle $p$ | N |
| $\mathbf{S}_u$ | Explicit momentum source vector | N/m³ |
| $\mathbf{U}$ | Fluid velocity field vector | m/s |
| $\mathbf{U}_c$ | Fluid velocity at cell center $c$ | m/s |
| $\mathbf{U}_f$ | Fluid velocity | m/s |
| $\mathbf{U}_{rel}$ | Relative velocity = $\mathbf{v}_p - \mathbf{U}_f$ | m/s |
| $\overline{\mathbf{U}}_f$ | Averaged (evaluated) fluid velocity | m/s |
| $\overline{\mathbf{U}}_p$ | Fluid velocity at particle location | m/s |
| $\mathbf{v}_p$ | Particle velocity vector | m/s |
| $\overline{\mathbf{v}}_c$ | Cell-averaged particle velocity | m/s |
| $\overline{\mathbf{v}}_p$ | Particle velocity in coupling calculations | m/s |
| $\mathbf{x}_i$ | Position vector of cell $i$ | m |
| $\mathbf{x}_p$ | Position vector of particle $p$ | m |

### Roman Symbols (Other)

| Symbol | Definition | Units |
|--------|-----------|-------|
| $\phi$ | Scalar field (general quantity) | variable |
| $w_i$ | Weight distribution for cell $i$ (Gaussian kernel) | - |
| $w_{p,c}$ | Weight of particle $p$ in cell $c$ (distribution method) | - |

### Subscripts and Superscripts

| Notation | Meaning |
|----------|---------|
| $(\cdot)_c$ | Property at or in cell $c$ |
| $(\cdot)_f$ | Property of fluid phase |
| $(\cdot)_p$ | Property of particle $p$ |
| $(\cdot)_{min}$ | Minimum value |
| $(\cdot)_{rel}$ | Relative or differential quantity |
| $(\cdot)_i$ | Component $i$ or cell/property $i$ |
| $\overline{(\cdot)}$ | Averaged or evaluated quantity |
| $(\cdot)^d$ | Dimensionless quantity |

### Dimensionless Numbers

| Symbol | Definition | Range/Application |
|--------|-----------|------------------|
| $Re$ | Reynolds number (particle) | - |
| $Re_p$ | Reynolds number (particle in lift) | Lift models |
| $Re_\omega$ | Rotation Reynolds number (particle spin) | Lift models |
| $Re_\Omega$ | Rotation Reynolds number (fluid vorticity) | Lift models |
| $Sr$ | Shear number | Lift models |
| $Rr$ | Rotation rate parameter | Lift models |

### Dimensional Specifications

#### Length
- m = meter
- 1/3 in exponent indicates cubic root

#### Mass
- kg = kilogram

#### Time
- s = second
- rad = radian

#### Compound Units
- m/s = meter per second (velocity)
- m²/s = square meter per second (diffusivity)
- m³ = cubic meter (volume)
- N = newton (force) = kg·m/s²
- Pa = pascal (pressure) = kg/(m·s²)
- Pa·s = pascal-second (dynamic viscosity)
- kg/m³ = kilogram per cubic meter (density)
- kg/(m³·s) = kilogram per cubic meter per second
- N/m³ = newton per cubic meter (force per volume)
- rad/s = radian per second (angular velocity)

---

## 10. References

1. Di Felice, R. (1994). The voidage function for fluid–particle interaction systems. International Journal of Multiphase Flow, 20(2), 153-159.
2. Gidaspow, D. (1994). Multiphase Flow and Fluidization: Continuum and Kinetic Theory Descriptions. Academic Press.
3. Beetstra, R., van der Hoef, M. A., & Kuipers, J. A. M. (2007). Drag force of intermediate Reynolds number flow past mono- and bidisperse arrays of spheres. AIChE Journal, 53(3), 489-501.
4. Loth, E. (2008). Lift of a spherical particle subject to vorticity and/or spin. AIAA Journal, 46(4), 801-809.
5. Rong, L. W., Dong, K. J., & Yu, A. B. (2013). Lattice-Boltzmann simulation of fluid flow through packed beds of uniform spheres. Chemical Engineering Science, 99, 44-58.
6. Shi, P., & Rzehak, R. (2019). Lift forces on solid spherical particles in unbounded flows. Chemical Engineering Science, 208, 115145.
7. Saffman, P. G. T. (1965). The lift on a small sphere in a slow shear flow. Journal of Fluid Mechanics, 22(2), 385-400.
8. Kianimoqadam, A., & Lapp, J. (2026). Gaussian integral method for void fraction. Particuology, 108, 125-142. [https://doi.org/10.1016/j.partic.2025.10.014](https://doi.org/10.1016/j.partic.2025.10.014)
