# Model Predictive Control of Spacecraft Angular Momentum using the Earth Magnetic Field

Agile missions require well controlled total angular momentum emerging from external torques such as aerodynamic drag or gravity gradient. Current solutions rely on an instantaneous reduction of the total angular momentum vector using magnetic torquers.

Due to the nature of magnetic desaturation with magnetic coils and the Earth magnetic field, the desaturation happens in the plane perpendicular to Earth’s magnetic field lines only. A model predictive approach could take the magnetic field lines of future orbit positions into account in order to find an optimized control command to minimize the magnitude of the total angular momentum error.

This repository contains a **buildable MATLAB + Simulink project** for magnetic-torquer (MTQ) momentum desaturation with:

- A **baseline instantaneous controller** (projection/pseudoinverse style).
- A **linear time-varying MPC** that uses **future B-field directions** along the orbit.
- A **Simulink model auto-builder** script (no manual block editing).
- A **full attitude + wheels + MTQ** simulation (LVLH reference) with **interactive 3D visualization**.
- **IGRF-14 magnetic field model** support (Aerospace Toolbox `igrfmagm`).

## Quick Start

**MATLAB (run this):**
- `run_matlab.m`

**Simulink (run this):**
- `run_simulink.m`

Both entry points automatically load parameters, run the selected model, and generate plots.

## How To Change Parameters

Default parameters live in:
- `MATLAB/Input/params/params_default.m`

You can either edit that file, or override in MATLAB before running:

```matlab
P = params_default();

% Model selection
P.plant.model = "full";        % "full" or "momentum"
P.controller  = "baseline";    % "baseline" or "mpc"
P.simulink.controller = "baseline"; % for Simulink runs

% Initial wheel momentum
P.x0 = [0.02; -0.01; 0.03];

% Orbit
P.orbit.alt_m = 550e3;
P.orbit.inc_deg = 97.5;

% Magnetic field model
P.env.bfield_model = "igrf";   % "igrf" or "dipole"
P.env.igrf_decimal_year = 2026.0;

% Run
run_matlab(P);      % MATLAB
run_simulink(P);    % Simulink
```

## Folder Structure

- `MATLAB/` (Input/Processing/Output)
- `Simulink/` (Input/Processing/Output)

**Detailed instructions:**
- `MATLAB/README.md`
- `Simulink/README.md`

## Requirements

- MATLAB (Optimization Toolbox for `quadprog`)
- Simulink (for the auto-built model)
- Aerospace Toolbox (for `igrfmagm`, `ecef2lla`, etc. when using IGRF)

If Aerospace Toolbox is not available, set `P.env.bfield_model = "dipole"`.
