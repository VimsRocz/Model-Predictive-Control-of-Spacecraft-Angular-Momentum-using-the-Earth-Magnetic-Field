# Model Predictive Control of Spacecraft Angular Momentum using the Earth Magnetic Field

Agile missions require well controlled total angular momentum emerging from external torques such as aerodynamic drag or gravity gradient. Current solutions rely on an instantaneous reduction of the total angular momentum vector using magnetic torquers.

Due to the nature of magnetic desaturation with magnetic coils and the Earth magnetic field, the desaturation happens in the plane perpendicular to Earth’s magnetic field lines only. A model predictive approach could take the magnetic field lines of future orbit positions into account in order to find an optimized control command to minimize the magnitude of the total angular momentum error.

This repository contains a **buildable MATLAB + Simulink project** for magnetic-torquer (MTQ) momentum desaturation with:

- A **baseline instantaneous controller** (projection/pseudoinverse style).
- A **linear time-varying MPC** that uses **future B-field directions** along the orbit.
- A **Simulink model auto-builder** script (no manual block editing).

## Quick start (MATLAB)

1. Open MATLAB at the repo root.
2. Run the simulation script:

	 - `main_run_sim.m`

This will generate `sim_out.mat` and open plots.

## Quick start (Simulink)

1. In MATLAB, run:

	 - `P = params_default();`
	 - `assignin('base','P',P);`
	 - `build_simulink_model();`
	 - `sim("mtq_desat_mpc_model");`

2. The signals `x_log` and `m_log` are saved to the workspace.

## Project structure

```
main_run_sim.m
build_simulink_model.m
params/
	params_default.m
models/
	skew.m
	orbit_propagator_simple.m
	ext_torques.m
controllers/
	baseline_mtq.m
	build_mpc_qp.m
	mpc_mtq_qp.m
plotting/
	plot_results.m
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

If you want a full attitude + wheel dynamics plant, or a real IGRF model, add your toolboxes and I will extend this project.
