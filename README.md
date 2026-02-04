# Model Predictive Control of Spacecraft Angular Momentum using the Earth Magnetic Field

Agile missions require well controlled total angular momentum emerging from external torques such as aerodynamic drag or gravity gradient. Current solutions rely on an instantaneous reduction of the total angular momentum vector using magnetic torquers.

Due to the nature of magnetic desaturation with magnetic coils and the Earth magnetic field, the desaturation happens in the plane perpendicular to Earth’s magnetic field lines only. A model predictive approach could take the magnetic field lines of future orbit positions into account in order to find an optimized control command to minimize the magnitude of the total angular momentum error.

This repository contains a **buildable MATLAB + Simulink project** for magnetic-torquer (MTQ) momentum desaturation with:

- A **baseline instantaneous controller** (projection/pseudoinverse style).
- A **linear time-varying MPC** that uses **future B-field directions** along the orbit.
- A **Simulink model auto-builder** script (no manual block editing).
- A **full attitude + wheels + MTQ** simulation (LVLH reference) with **interactive 3D visualization**.
- **IGRF-14 magnetic field model** support (Aerospace Toolbox `igrfmagm`).

## Quick start (MATLAB)

1. Open MATLAB at the repo root.
2. Run the simulation script:

	 - `main_run_sim.m`

This will generate `sim_out.mat` and open plots.

### Full attitude + wheels plant (MATLAB)

- `main_run_full_sim("baseline")`
- `main_run_full_sim("mpc")`

Outputs are saved to `outputs/matlab_full/` and include a 3D scene that visualizes the **“MTQ torque is always ⟂B”** constraint using a plane at the satellite position.

## Quick start (Simulink)

1. In MATLAB, run:

	 - `P = params_default();`
	 - `assignin('base','P',P);`
	 - `build_simulink_model();`
	 - `sim("mtq_desat_mpc_model");`

2. The signals `x_log` and `m_log` are saved to the workspace.

> Note: if `P.env.bfield_model = "igrf"`, `run_simulink.m` automatically uses an **Interpreted MATLAB Function** variant of the model so it can call `igrfmagm`.

### Full plant in Simulink

- `run_simulink_full(P)`

This builds `mtq_full_model.slx` (with an interpreted MATLAB Function block that calls `sim_full_step`) and writes results to `outputs/simulink_full/`.

## Project structure

```
main_run_sim.m
main_run_full_sim.m
build_simulink_model.m
build_simulink_model_imf.m
interactive_mtq_explorer.m
params/
	params_default.m
models/
	skew.m
	earth_dipole_field.m
	earth_igrf_field.m
	earth_mag_field_eci.m
	orbit_propagator_simple.m
	orbit_state_simple.m
	orbit_state_rv_simple.m
	ext_torques.m
	ext_torques_body.m
	lvlh_reference.m
controllers/
	baseline_mtq.m
	build_mpc_qp.m
	mpc_mtq_qp.m
plotting/
	plot_results.m
	plot_results_full.m
```

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
