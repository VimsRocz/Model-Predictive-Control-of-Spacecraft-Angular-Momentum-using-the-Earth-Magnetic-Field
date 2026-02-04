# MATLAB

## INPUT
- Default parameters: `MATLAB/Input/params/params_default.m`
- To change parameters, edit the file above, or use a custom struct in MATLAB:

```matlab
P = params_default();
P.plant.model = "full";        % "full" or "momentum"
P.controller  = "baseline";    % "baseline" or "mpc"
P.x0 = [0.02; -0.01; 0.03];      % initial wheel momentum
P.orbit.alt_m = 550e3;
P.orbit.inc_deg = 97.5;
P.env.bfield_model = "igrf";   % "igrf" or "dipole"
P.env.igrf_decimal_year = 2026.0;

% 3D visualization controls
P.viz.show_sun = true;          % show Sun in 3D scene
P.viz.animate_all_steps = true; % full rotation animation
P.viz.animate_stride = 1;       % use >1 to speed up animation
P.viz.auto_animate_full = true; % auto-play after run
P.viz.auto_animate_momentum = true; % auto-play (momentum-only) after run
P.viz.interactive = true;       % rotate/pan/zoom tools
P.viz.use_lighting = true;      % 3D lighting for Earth/Sun
P.viz.earth_texture = true;     % textured Earth (built-in topo) when available

run_matlab(P);
```

## PROCESSING
- Main entry point (run this): `run_matlab.m` (repo root)
- Internal drivers (do not run directly):
  - `MATLAB/Processing/main_run_sim.m` (momentum-only)
  - `MATLAB/Processing/main_run_full_sim.m` (full attitude + wheels)
- Core models/controllers live under:
  - `MATLAB/Processing/models/`
  - `MATLAB/Processing/controllers/`
  - `MATLAB/Processing/plotting/`

## OUTPUT
- Momentum-only outputs:
  - `MATLAB/Output/momentum/sim_out_matlab.mat`
  - `MATLAB/Output/momentum/report_timeseries_*.png`
  - `MATLAB/Output/momentum/scene_earth_*.png`
- Full model outputs:
  - `MATLAB/Output/full/sim_out_full.mat`
  - `MATLAB/Output/full/report_full_timeseries_*.png`
  - `MATLAB/Output/full/scene_full_earth_*.png`
- Benchmarks (optional):
  - `MATLAB/Output/analysis/benchmark_summary.*`
