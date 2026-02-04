function tau = sim_tau_ext(t)
%SIM_TAU_EXT Wrapper for Simulink Interpreted MATLAB Function block.

persistent P
if isempty(P)
    P = evalin('base','P');
end

k = floor(t/P.Ts) + 1;
tau = ext_torques(k, P);
end
