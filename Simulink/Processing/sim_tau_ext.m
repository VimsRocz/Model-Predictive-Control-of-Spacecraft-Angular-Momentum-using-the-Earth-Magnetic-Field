function tau = sim_tau_ext(t)
%SIM_TAU_EXT Wrapper for Simulink Interpreted MATLAB Function block.

persistent P last_run_id
[baseP, base_id] = load_base_P();
if isempty(P) || isempty(last_run_id) || ~strcmp(last_run_id, base_id)
    P = baseP;
    last_run_id = base_id;
end

k = floor(t/P.Ts) + 1;
tau = ext_torques(k, P);
end

function [Pbase, run_id] = load_base_P()
Pbase = evalin('base','P');
run_id = "";
if isfield(Pbase,'simulink') && isfield(Pbase.simulink,'run_id')
    run_id = string(Pbase.simulink.run_id);
end
run_id = char(run_id);
end
