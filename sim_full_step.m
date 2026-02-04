function [r_eci, B_eci, B_body, h_w, m, tau_ext, tau_mtq, tau_rw, q_ib, q_ib_des, w_body, w_des_body] = sim_full_step(t)
%SIM_FULL_STEP One-step full-plant update for Simulink (Interpreted MF).
%
% Uses persistent state so Simulink can call this as a single block.
% Outputs are the *current* state at time t (not the next state).

persistent P state last_k last_t

if isempty(P)
    P = evalin('base','P');
end

if isempty(last_t) || t < last_t
    [state, last_k] = init_state(P);
end
last_t = t;

k = floor(t / P.Ts) + 1;
if k < 1
    k = 1;
end

% Advance state to the requested step if needed.
if isempty(last_k)
    [state, last_k] = init_state(P);
end
if k > last_k
    for kk = last_k:(k-1)
        state = advance_one_step(state, kk, P);
    end
    last_k = k;
end

% Compute outputs at step k using current state.
[r_eci, B_eci, B_body, h_w, m, tau_ext, tau_mtq, tau_rw, q_ib, q_ib_des, w_body, w_des_body] = ...
    compute_outputs(state, k, P);
end

% ---------------- helpers ----------------
function [state, k] = init_state(P)
[r0, v0] = orbit_state_rv_simple(1, P);
[q0_des, w0_des] = lvlh_reference(r0, v0);

state.q = quatnormalize(q0_des);
state.w = w0_des(:);

if isfield(P,'x0') && numel(P.x0) == 3
    state.h = P.x0(:);
else
    state.h = [0.02; -0.01; 0.03];
end

state.m_prev = zeros(3,1);
state.U_prev = [];
state.m_hold = zeros(3,1);

state.last_solve = 0;
k = 1;
end

function state = advance_one_step(state, k, P)
t_s = (k-1) * P.Ts;

[r_eci, v_eci] = orbit_state_rv_simple(k, P);
[q_des, w_des] = lvlh_reference(r_eci, v_eci);

q = state.q;
w = state.w;
h = state.h;

C_bi = quat2dcm(q);
Bk_eci = earth_mag_field_eci(r_eci, P, t_s);
Bk_body = C_bi * Bk_eci;

tau_ext = ext_torques_body(t_s, r_eci, v_eci, q, P);

% Desaturation command
mode = "baseline";
if isfield(P,'simulink') && isfield(P.simulink,'controller')
    mode = string(P.simulink.controller);
end
mode = lower(mode);

if mode == "baseline"
    mk = baseline_mtq(h, Bk_body, P);
else
    [mk, state] = mpc_command_for_step(state, k, P);
end

tau_mtq = cross(mk, Bk_body);

% Integration (internal substeps)
nSub = 10;
if isfield(P,'full') && isfield(P.full,'int_substeps')
    nSub = max(1, round(P.full.int_substeps));
end
dt_int = P.Ts / nSub;

I = P.sc.I_kgm2;

Kp = 0.05; Kd = 0.2;
if isfield(P,'adcs') && isfield(P.adcs,'Kp'); Kp = double(P.adcs.Kp); end
if isfield(P,'adcs') && isfield(P.adcs,'Kd'); Kd = double(P.adcs.Kd); end

tau_rw_max = inf(3,1);
if isfield(P,'rw') && isfield(P.rw,'tau_max_Nm'); tau_rw_max = P.rw.tau_max_Nm(:); end

for sub = 1:nSub
    q_err = quatmultiply(quatconj(q_des), q);
    if q_err(1) < 0
        q_err = -q_err;
    end
    e = q_err(2:4).';
    w_err = w - w_des(:);
    tau_rw = -Kp * e - Kd * w_err;
    tau_rw = min(max(tau_rw, -tau_rw_max), tau_rw_max);

    h_dot = -tau_rw;
    w_dot = I \ (tau_ext + tau_mtq + tau_rw - cross(w, I*w + h));
    q_dot = 0.5 * quatmultiply(q, [0 w.']);

    h = h + dt_int * h_dot;
    w = w + dt_int * w_dot;
    q = quatnormalize(q + dt_int * q_dot);
end

state.q = q;
state.w = w;
state.h = h;
state.m_prev = mk;
end

function [r_eci, B_eci, B_body, h_w, m, tau_ext, tau_mtq, tau_rw, q_ib, q_ib_des, w_body, w_des_body] = ...
    compute_outputs(state, k, P)
t_s = (k-1) * P.Ts;

[r_eci, v_eci] = orbit_state_rv_simple(k, P);
[q_des, w_des] = lvlh_reference(r_eci, v_eci);

q = state.q;
w = state.w;
h = state.h;

C_bi = quat2dcm(q);

B_eci = earth_mag_field_eci(r_eci, P, t_s);
B_body = C_bi * B_eci;

tau_ext = ext_torques_body(t_s, r_eci, v_eci, q, P);

% Desaturation command (current)
mode = "baseline";
if isfield(P,'simulink') && isfield(P.simulink,'controller')
    mode = string(P.simulink.controller);
end
mode = lower(mode);

if mode == "baseline"
    m = baseline_mtq(h, B_body, P);
else
    [m, ~] = mpc_command_for_step(state, k, P);
end

tau_mtq = cross(m, B_body);

% PD wheel torque (current)
Kp = 0.05; Kd = 0.2;
if isfield(P,'adcs') && isfield(P.adcs,'Kp'); Kp = double(P.adcs.Kp); end
if isfield(P,'adcs') && isfield(P.adcs,'Kd'); Kd = double(P.adcs.Kd); end

tau_rw_max = inf(3,1);
if isfield(P,'rw') && isfield(P.rw,'tau_max_Nm'); tau_rw_max = P.rw.tau_max_Nm(:); end

q_err = quatmultiply(quatconj(q_des), q);
if q_err(1) < 0
    q_err = -q_err;
end
e = q_err(2:4).';
w_err = w - w_des(:);
tau_rw = -Kp * e - Kd * w_err;
tau_rw = min(max(tau_rw, -tau_rw_max), tau_rw_max);

% Outputs
r_eci = r_eci(:);
B_eci = B_eci(:);
B_body = B_body(:);
h_w = h(:);
m = m(:);
tau_ext = tau_ext(:);
tau_mtq = tau_mtq(:);
tau_rw = tau_rw(:);
q_ib = q(:).';
q_ib_des = q_des(:).';
w_body = w(:);
w_des_body = w_des(:);
end

function [m_cmd, state] = mpc_command_for_step(state, k, P)
% Shared MPC handler (uses same logic as MATLAB sim).

update_every = 1;
warm_start = false;
if isfield(P,'mpc') && isfield(P.mpc,'update_every')
    update_every = max(1, round(P.mpc.update_every));
end
if isfield(P,'mpc') && isfield(P.mpc,'warm_start')
    warm_start = logical(P.mpc.warm_start);
end

doSolve = (mod(k-1, update_every) == 0) || (k == 1);
if doSolve
    U0 = [];
    if warm_start && ~isempty(state.U_prev)
        U0 = [state.U_prev(4:end); state.U_prev(end-2:end)];
    end

    Bprov = @(kk) pred_B_body(kk, P);
    tauprov = @(kk) pred_tau_ext_body(kk, P);
    [m_cmd, U_opt] = mpc_mtq_qp(state.h, k, P, Bprov, tauprov, state.m_prev, U0);
    state.U_prev = U_opt;
    state.m_hold = m_cmd;
else
    m_cmd = state.m_hold;
end
end

function Bk = pred_B_body(kk, P)
kk = max(1, kk);
t_s = (kk-1) * P.Ts;
[r, v] = orbit_state_rv_simple(kk, P);
[~, ~, C_bi_des] = lvlh_reference(r, v);
Bk_eci = earth_mag_field_eci(r, P, t_s);
Bk = C_bi_des * Bk_eci;
end

function tauk = pred_tau_ext_body(kk, P)
kk = max(1, kk);
t_s = (kk-1) * P.Ts;
[r, v] = orbit_state_rv_simple(kk, P);
[q_des, ~] = lvlh_reference(r, v);
tauk = ext_torques_body(t_s, r, v, q_des, P);
end

