function S = simulate_mtq_full(P, desatMode)
%SIMULATE_MTQ_FULL Full attitude + wheels + MTQ desaturation simulation.
%
% States (body frame unless noted):
%   q_ib  : quaternion inertial->body (scalar-first)
%   w     : body angular rate [rad/s]
%   h_w   : reaction-wheel momentum [N*m*s]
%
% Dynamics:
%   q_dot = 0.5 * q ⊗ [0; w]
%   I*w_dot = tau_ext + tau_mtq + tau_rw - w×(I*w + h_w)
%   h_w_dot = -tau_rw   (internal wheel torque exchanges momentum with body)
%
% MTQ torque:
%   tau_mtq = m × B_body
%
% Desaturation controller:
%   baseline: tau_des = -KH*h_w, m = (B×tau_des)/||B||^2 with saturation
%   mpc:     uses mpc_mtq_qp() on a wheel-momentum proxy model (same as momentum-only)
%
% Attitude reference:
%   LVLH (nadir-pointing) using lvlh_reference().

if nargin < 2 || strlength(desatMode) == 0
    desatMode = "baseline";
end
desatMode = lower(string(desatMode));

% Pre-allocate logs
N = P.N;
t = (0:N-1) * P.Ts;

r_eci = zeros(N,3);
v_eci = zeros(N,3);

q_ib = zeros(N,4);
w_body = zeros(N,3);
h_w = zeros(N,3);

q_ib_des = zeros(N,4);
w_des_body = zeros(N,3);

B_eci = zeros(N,3);
B_body = zeros(N,3);

m_cmd = zeros(N,3);

tau_ext_body = zeros(N,3);
tau_mtq_body = zeros(N,3);
tau_rw_body = zeros(N,3);

% ADCS gains / limits (defaults if not provided)
Kp = 0.05;
Kd = 0.2;
if isfield(P,'adcs') && isfield(P.adcs,'Kp'); Kp = double(P.adcs.Kp); end
if isfield(P,'adcs') && isfield(P.adcs,'Kd'); Kd = double(P.adcs.Kd); end

tau_rw_max = inf(3,1);
if isfield(P,'rw') && isfield(P.rw,'tau_max_Nm'); tau_rw_max = P.rw.tau_max_Nm(:); end

I = P.sc.I_kgm2;

% Internal integration sub-steps (improves stability vs 1 s Euler updates)
nSub = 10;
if isfield(P,'full') && isfield(P.full,'int_substeps')
    nSub = max(1, round(P.full.int_substeps));
end
dt_int = P.Ts / nSub;

% Initial orbit, reference attitude, and state
[r0, v0] = orbit_state_rv_simple(1, P);
[q0_des, w0_des] = lvlh_reference(r0, v0);

q_ib(1,:) = quatnormalize(q0_des);
w_body(1,:) = w0_des(:).';

if isfield(P,'x0') && numel(P.x0) == 3
    h_w(1,:) = P.x0(:).';
else
    h_w(1,:) = [0.02 -0.01 0.03];
end

% MPC plumbing
m_prev = zeros(3,1);
U_prev = [];
m_hold = zeros(3,1);

update_every = 1;
warm_start = false;
if isfield(P,'mpc') && isfield(P.mpc,'update_every'); update_every = max(1, round(P.mpc.update_every)); end
if isfield(P,'mpc') && isfield(P.mpc,'warm_start'); warm_start = logical(P.mpc.warm_start); end

% Prediction providers (use desired LVLH attitude to keep model "one-frame")
B_provider   = @(kk) pred_B_body(kk);
tau_provider = @(kk) pred_tau_ext_body(kk);

for k = 1:N
    tk = t(k);

    % Orbit state (ECI)
    [rk, vk] = orbit_state_rv_simple(k, P);
    r_eci(k,:) = rk.';
    v_eci(k,:) = vk.';

    % Desired LVLH reference
    [qk_des, wk_des] = lvlh_reference(rk, vk);
    q_ib_des(k,:) = quatnormalize(qk_des);
    w_des_body(k,:) = wk_des(:).';

    % Current attitude DCM (inertial->body)
    qk = q_ib(k,:);
    C_bi = quat2dcm(qk);

    % Magnetic field (ECI -> body)
    Bk_eci = earth_mag_field_eci(rk, P, tk);
    B_eci(k,:) = Bk_eci.';
    Bk_body = C_bi * Bk_eci;
    B_body(k,:) = Bk_body.';

    % External torques in body frame (uses q)
    tau_ext_k = ext_torques_body(tk, rk, vk, qk, P);
    tau_ext_body(k,:) = tau_ext_k.';

    % Attitude control (wheels): PD on quaternion error (computed at sample for logging).
    q_err0 = quatmultiply(quatconj(q_ib_des(k,:)), qk); % q_des^{-1} ⊗ q
    if q_err0(1) < 0
        q_err0 = -q_err0;
    end
    e0 = q_err0(2:4).'; % 3x1
    w_err0 = (w_body(k,:).' - wk_des(:));
    tau_rw_k = -Kp * e0 - Kd * w_err0;
    tau_rw_k = min(max(tau_rw_k, -tau_rw_max), tau_rw_max);
    tau_rw_body(k,:) = tau_rw_k.';

    % Desaturation controller (MTQs) acts on wheel momentum h_w
    xk = h_w(k,:).';
    if desatMode == "baseline"
        mk = baseline_mtq(xk, Bk_body, P);
    elseif desatMode == "mpc"
        doSolve = (mod(k-1, update_every) == 0) || (k == 1);
        if doSolve
            U0 = [];
            if warm_start && ~isempty(U_prev)
                U0 = [U_prev(4:end); U_prev(end-2:end)];
            end
            [mk, U_opt] = mpc_mtq_qp(xk, k, P, B_provider, tau_provider, m_prev, U0);
            U_prev = U_opt;
            m_hold = mk;
        else
            mk = m_hold;
        end
    else
        error("Unknown desatMode '%s'. Use 'baseline' or 'mpc'.", desatMode);
    end

    m_cmd(k,:) = mk.';
    tau_mtq_k = cross(mk, Bk_body);
    tau_mtq_body(k,:) = tau_mtq_k.';

    % Integrate dynamics
    if k < N
        wk = w_body(k,:).';
        hk = h_w(k,:).';
        qint = qk;

        % Hold MTQ command constant over this sample.
        for sub = 1:nSub
            % Recompute PD wheel torque at the internal integration rate (improves stability).
            q_err = quatmultiply(quatconj(q_ib_des(k,:)), qint);
            if q_err(1) < 0
                q_err = -q_err;
            end
            e = q_err(2:4).';
            w_err = wk - wk_des(:);
            tau_rw = -Kp * e - Kd * w_err;
            tau_rw = min(max(tau_rw, -tau_rw_max), tau_rw_max);

            % Wheel momentum exchange
            h_dot = -tau_rw;

            % Body rotational dynamics
            w_dot = I \ (tau_ext_k + tau_mtq_k + tau_rw - cross(wk, I*wk + hk));

            % Quaternion kinematics
            q_dot = 0.5 * quatmultiply(qint, [0 wk.']);

            hk = hk + dt_int * h_dot;
            wk = wk + dt_int * w_dot;
            qint = quatnormalize(qint + dt_int * q_dot);
        end

        h_w(k+1,:) = hk.';
        w_body(k+1,:) = wk.';
        q_ib(k+1,:) = qint;
    end

    m_prev = mk;
end

S = struct();
S.P = P;
S.desatMode = desatMode;
S.t = t;
S.r_eci = r_eci;
S.v_eci = v_eci;
S.q_ib = q_ib;
S.w_body = w_body;
S.h_w = h_w;
S.q_ib_des = q_ib_des;
S.w_des_body = w_des_body;
S.B_eci = B_eci;
S.B_body = B_body;
S.m = m_cmd;
S.tau_ext_body = tau_ext_body;
S.tau_mtq_body = tau_mtq_body;
S.tau_rw_body = tau_rw_body;

% Helpful derived (in body)
S.tau_total_body = tau_ext_body + tau_mtq_body + tau_rw_body;
S.tau_mtq_dotB = sum(tau_mtq_body .* (B_body ./ max(vecnorm(B_body,2,2), 1e-30)), 2);

% ---------------- nested helper providers ----------------
function Bk = pred_B_body(kk)
    kk = min(max(kk, 1), N);
    tt = (kk-1) * P.Ts;
    [rr, vv] = orbit_state_rv_simple(kk, P);
    [q_des, ~, C_bi_des] = lvlh_reference(rr, vv);
    %#ok<NASGU> q_des
    Bk_eci2 = earth_mag_field_eci(rr, P, tt);
    Bk = C_bi_des * Bk_eci2;
end

function tauk = pred_tau_ext_body(kk)
    kk = min(max(kk, 1), N);
    tt = (kk-1) * P.Ts;
    [rr, vv] = orbit_state_rv_simple(kk, P);
    [q_des, ~] = lvlh_reference(rr, vv);
    tauk = ext_torques_body(tt, rr, vv, q_des, P);
end
end
