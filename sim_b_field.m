function B = sim_b_field(t)
%SIM_B_FIELD Wrapper for Simulink Interpreted MATLAB Function block.

persistent P
if isempty(P)
    P = evalin('base','P');
end

k = floor(t/P.Ts) + 1;
B = orbit_propagator_simple(k, P);
end
