# Laminar Flow Between Parallel Plates — CFD Validation and Thermal Extension

A computational fluid dynamics study of low-Reynolds-number flow through a two-dimensional parallel-plate channel. The project combines an analytical Poiseuille-flow solution, ANSYS Fluent simulations, a three-level mesh comparison, outlet-boundary-condition testing, and a non-isothermal extension.

## Problem statement

The objective is to model steady, incompressible, laminar flow through a straight channel bounded by two stationary, parallel walls. The numerical velocity profiles are compared with the fully developed analytical solution to assess how mesh resolution and outlet boundary-condition choice affect the predicted flow field. A separate thermal case examines the temperature field when the two walls are maintained at different temperatures.

## Geometry and physical parameters

| Parameter | Value |
|---|---:|
| Channel height, *H* | 0.020 m |
| Channel length, *L* | 0.100 m |
| Reynolds number, *Re* | 12 |
| Density, ρ | 1000 kg m⁻³ |
| Dynamic viscosity, μ | 1.02 × 10⁻³ Pa s |
| Kinematic viscosity, ν = μ/ρ | 1.02 × 10⁻⁶ m² s⁻¹ |
| Hydraulic diameter, *D*ₕ = 2*H* | 0.040 m |
| Mean inlet/reference velocity | 3.06 × 10⁻⁴ m s⁻¹ |
| Thermal diffusivity, α | 1.43 × 10⁻⁷ m² s⁻¹ |
| Specific heat capacity, *c*ₚ | 4185 J kg⁻¹ K⁻¹ |
| Derived thermal conductivity, *k* = αρ*c*ₚ | 0.598455 W m⁻¹ K⁻¹ |

The thermal extension uses an inlet temperature of 300 K, a lower-wall temperature of 50 K, and an upper-wall temperature of 350 K.

## Analytical reference

For fully developed, incompressible, laminar flow between stationary parallel plates, with the transverse coordinate *y* measured from the channel centreline, the streamwise velocity is

\[
 u(y)=\frac{\Delta p}{2\mu L}\left(\frac{H^2}{4}-y^2\right),
 \qquad -\frac{H}{2}\le y\le\frac{H}{2}.
\]

The corresponding mean velocity is

\[
 U_{\mathrm{mean}}=\frac{\Delta p H^2}{12\mu L},
 \qquad u_{\max}=\frac{3}{2}U_{\mathrm{mean}}.
\]

The MATLAB post-processing converts between this centreline coordinate and the wall-based coordinate used in the Fluent exports. With *Y* = 0 at the lower wall, the equivalent profile is

\[
 u(Y)=6U_{\mathrm{mean}}\frac{Y}{H}\left(1-\frac{Y}{H}\right).
\]

For the parameters above, the ideal fully developed pressure drop is approximately 9.3636 × 10⁻⁴ Pa.

## Computational approach

1. A two-dimensional channel model is prepared and solved in ANSYS Workbench / Fluent.
2. Three structured mesh resolutions are used to examine spatial-resolution effects.
3. Outlet predictions are compared using **Pressure Outlet** and **Outflow** conditions.
4. Outlet X-velocity profiles are sampled at Fluent cell-centre coordinates and compared against the analytical profile.
5. Mean absolute error, maximum absolute error, and normalised discrete L₂ error are calculated.
6. First- and second-order upwind schemes are compared on the coarse mesh.
7. A thermal case with unequal wall temperatures is post-processed separately.

### Mesh set

| Mesh | Nominal element size | Cells | Minimum orthogonal quality | Maximum aspect ratio |
|---|---:|---:|---:|---:|
| Coarse | 0.00050 m | 8,000 | 0.999691 | 1.44151 |
| Medium | 0.00025 m | 32,000 | 0.998463 | 1.46573 |
| Fine | 0.000125 m | 128,000 | 0.998801 | 1.45933 |

All three meshes have high reported orthogonal quality and low aspect ratios. The available error values are not monotonic with refinement, so the results are presented as a mesh-comparison study rather than a claim of formal grid independence.

## Selected results

The medium mesh has the lowest normalised L₂ velocity error among the three tested meshes for both outlet conditions:

| Outlet condition | Coarse | Medium | Fine |
|---|---:|---:|---:|
| Pressure Outlet | 0.13334% | **0.03677%** | 0.22076% |
| Outflow | 0.34530% | **0.09125%** | 0.14232% |

For the coarse-mesh scheme comparison, the normalised L₂ errors are 0.13369% for first-order upwind and 0.13334% for second-order upwind. The difference is small for this sampled profile.

The thermal case includes outlet-temperature data, static-temperature contours, and supporting convergence and flow-field images. The repository retains the exported data used by the MATLAB analysis so the reported comparisons can be regenerated.

## Visual results

### Mesh resolutions

![Medium mesh](results/mesh/medium_mesh.png)

### Axial velocity field — medium mesh, Pressure Outlet

![Axial velocity contour](results/contours/X_Velocity_medium_P.png)

### Temperature field — thermal case

![Temperature contour](results/contours/static_Temp_medium_P_thermal.png)

## Repository structure

~~~text
.
├── README.md
├── ansys/
│   └── parallel_plate_channel_flow.wbpz
├── matlab/
│   ├── analyse_parallel_plate_flow.m
│   └── figures_matlab/          # MATLAB .fig and 300 dpi .png outputs
├── data/                        # Fluent profile exports and CSV summaries
└── results/
    ├── mesh/                    # Coarse, medium and fine mesh images
    ├── contours/                # Velocity, pressure and temperature contours
    ├── residuals/               # Fluent residual-history images
    └── animations/              # Flow visualisation video
~~~

## Reproducing the MATLAB post-processing

1. Open MATLAB.
2. Set the **Current Folder** to this repository's `data` directory.
3. Run the following command in the MATLAB Command Window:

~~~matlab
run(fullfile('..', 'matlab', 'analyse_parallel_plate_flow.m'))
~~~

The script reads the Fluent-exported profile files from the current folder. It regenerates the analytical reference values, pointwise error tables, comparison summaries, mesh statistics, and available plots. New figures are saved to `data/figures_matlab/` as PNG and MATLAB FIG files.

The thermal outlet-temperature and centreline-pressure plots are generated when their corresponding input files are present. The Fluent residual images are retained as solver-output evidence and are not generated by the MATLAB script.

## Scope and interpretation

- The velocity comparison uses the exported cell-centre samples, not an independently reconstructed continuous CFD field.
- Error values are specific to the supplied mesh, boundary-condition, and solver data.
- The non-monotonic error trend means that the finest mesh does not produce the smallest error for every case. No formal Richardson extrapolation or Grid Convergence Index is claimed.
- The ANSYS Workbench archive is provided as a project package; its internal solver states and compatibility depend on the ANSYS release used to open it.

## Software

- ANSYS Workbench / Fluent — geometry, meshing, flow and thermal simulations
- MATLAB — analytical reference, post-processing, error calculations and plots
