# Model Predictive Control of Spacecraft Angular Momentum using the Earth Magnetic Field

Agile missions require well controlled total angular momentum emerging from external torques such as aerodynamic drag or gravity gradient. Current solutions rely on an instantaneous reduction of the total angular momentum vector using magnetic torquers.

Due to the nature of magnetic desaturation with magnetic coils and the Earth magnetic field, the desaturation happens in the plane perpendicular to Earth’s magnetic field lines only. A model predictive approach could take the magnetic field lines of future orbit positions into account in order to find an optimized control command to minimize the magnitude of the total angular momentum error.

This repository contains a **buildable MATLAB + Simulink project** for magnetic-torquer (MTQ) momentum desaturation with:

- A **baseline instantaneous controller** (projection/pseudoinverse style).
- A **linear time-varying MPC** that uses **future B-field directions** along the orbit.
- A **Simulink model auto-builder** script (no manual block editing).
- A **full attitude + wheels + MTQ** simulation (LVLH reference) with **interactive 3D visualization**.
- **IGRF-14 magnetic field model** support (Aerospace Toolbox `igrfmagm`).

## **Input**

**Run files (only these two are entry points):**
- `run_matlab.m`
- `run_simulink.m`

**Parameter file:**
- `MATLAB/Input/params/params_default.m`

**Folder structure (requested):**
- `MATLAB/Input`, `MATLAB/Processing`, `MATLAB/Output`
- `Simulink/Input`, `Simulink/Processing`, `Simulink/Output`

**Common parameters (override via UI or by passing a custom `P`):**
- `P.plant.model = "full" | "momentum"`
- `P.controller = "baseline" | "mpc"`
- `P.x0 = [h_x; h_y; h_z]` (initial wheel momentum / x0)
- `P.orbit.alt_m`, `P.orbit.inc_deg`
- `P.env.bfield_model = "igrf" | "dipole"`
- `P.env.igrf_decimal_year`

**Interactive UI (optional):**
- `interactive_mtq_explorer.m`

## **Process**

**MATLAB pipeline (called by `run_matlab.m`):**
- `simulate_mtq_full.m` (full model) or `simulate_mtq.m` (momentum-only)
- `controllers/baseline_mtq.m`
- `controllers/mpc_mtq_qp.m`
- `models/earth_mag_field_eci.m` → `models/earth_igrf_field.m` or `models/earth_dipole_field.m`
- `models/ext_torques_body.m` (full) or `models/ext_torques.m` (momentum-only)
- Internal drivers: `main_run_full_sim.m`, `main_run_sim.m` (do not run directly)
- Logical grouping: see `MATLAB/Processing/README.md`

**Simulink pipeline (called by `run_simulink.m`):**
- Full model: `build_simulink_full_model_imf.m` → `mtq_full_model.slx`
- Full model stepper: `sim_full_step.m`
- Momentum-only model: `build_simulink_model.m` or `build_simulink_model_imf.m`
- IMF helpers: `sim_b_field.m`, `sim_controller.m`, `sim_tau_ext.m`
- Logical grouping: see `Simulink/Processing/README.md`

## **Output**

**MATLAB outputs (full model):**
- `MATLAB/Output/full/sim_out_full.mat`
- `MATLAB/Output/full/report_full_timeseries_*.png`
- `MATLAB/Output/full/scene_full_earth_*.png`

**MATLAB outputs (momentum-only):**
- `MATLAB/Output/momentum/sim_out_matlab.mat`
- `MATLAB/Output/momentum/report_timeseries_*.png`
- `MATLAB/Output/momentum/scene_earth_*.png`

**MATLAB outputs (benchmarks):**
- `MATLAB/Output/analysis/benchmark_summary.*`

**Simulink outputs (full model):**
- `Simulink/Output/full/sim_out_full_simulink.mat`
- `Simulink/Output/full/report_full_timeseries_Simulink.png`
- `Simulink/Output/full/scene_full_earth_Simulink.png`

**Simulink outputs (momentum-only):**
- `Simulink/Output/momentum/sim_out_simulink.mat`
- `Simulink/Output/momentum/report_timeseries_Simulink.png`
- `Simulink/Output/momentum/scene_earth_Simulink.png`

**MATLAB vs Simulink consistency:**
- `run_matlab.m` and `run_simulink.m` now use the same controller (`P.controller` / `P.simulink.controller`) and the same plant (`P.plant.model`). Both print a max-difference summary if the other output exists.

## Folders

- `MATLAB/` (Input/Processing/Output)
- `Simulink/` (Input/Processing/Output)
- `models/` dynamics, orbit, and environment models
- `controllers/` baseline + MPC controllers
- `plotting/` 2D reports and 3D scenes/animations

## Controller notes

**Baseline:**

```
m = (B × τ_des) / ||B||^2,  τ_des = -K_H x
```

**MPC:** solves a box-constrained QP over a horizon with time-varying `B_k`.

## Benchmark ideas

- Horizon length sweep (10 / 30 / 60 / 120)
- MTQ authority sweep (`m_max` low vs high)
- Orbit inclination changes (affects `B` geometry)
- Bias + periodic disturbance torque
- Monte Carlo on inertia and disturbance uncertainty

## Requirements

- MATLAB (Optimization Toolbox for `quadprog`)
- Simulink (for the auto-built model)
- Aerospace Toolbox (for `igrfmagm`, `ecef2lla`, etc. when using IGRF)

If Aerospace Toolbox is not available, set `P.env.bfield_model = "dipole"`.
