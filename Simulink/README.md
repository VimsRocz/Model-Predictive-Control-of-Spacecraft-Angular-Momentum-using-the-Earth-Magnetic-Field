# Simulink

## INPUT
- Parameters are shared with MATLAB:
  - `MATLAB/Input/params/params_default.m`
- Change parameters in the file above, or pass a custom struct:

```matlab
P = params_default();
P.plant.model = "full";              % "full" or "momentum"
P.simulink.controller = "baseline"; % "baseline" or "mpc"
P.orbit.alt_m = 550e3;
P.orbit.inc_deg = 97.5;
P.env.bfield_model = "igrf";         % "igrf" or "dipole"
P.env.igrf_decimal_year = 2026.0;

% 3D visualization controls
P.viz.show_sun = false;         % show Sun in 3D scene (optional)
P.viz.animate_all_steps = true; % full rotation animation
P.viz.animate_stride = 1;       % use >1 to speed up animation
P.viz.auto_animate_full = true; % auto-play after run
P.viz.auto_animate_momentum = true; % auto-play (momentum-only) after run
P.viz.interactive = true;       % rotate/pan/zoom tools
P.viz.use_lighting = true;      % 3D lighting for Earth
P.viz.earth_texture = true;     % textured Earth (built-in topo) when available

run_simulink(P);
```

## PROCESSING
- Main entry point (run this): `run_simulink.m` (repo root)
- Simulink builders + IMF wrappers live under:
  - `Simulink/Processing/build_simulink_model.m`
  - `Simulink/Processing/build_simulink_model_imf.m`
  - `Simulink/Processing/build_simulink_full_model_imf.m`
  - `Simulink/Processing/sim_b_field.m`, `sim_controller.m`, `sim_tau_ext.m`, `sim_full_step.m`

## OUTPUT
- Momentum-only outputs:
  - `Simulink/Output/momentum/sim_out_simulink.mat`
  - `Simulink/Output/momentum/report_timeseries_Simulink.png`
  - `Simulink/Output/momentum/scene_earth_Simulink.png`
- Full model outputs:
  - `Simulink/Output/full/sim_out_full_simulink.mat`
  - `Simulink/Output/full/report_full_timeseries_Simulink.png`
  - `Simulink/Output/full/scene_full_earth_Simulink.png`
