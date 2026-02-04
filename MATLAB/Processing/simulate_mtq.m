function S = simulate_mtq(P, ctrlMode)
%SIMULATE_MTQ Momentum-only MTQ desaturation simulation (baseline or MPC).
%
% State (momentum error):
%   x_{k+1} = x_k + Ts*(tau_ext + tau_mtq)
%   tau_mtq = cross(m, B) = -skew(B)*m
%
% Notes:
%   - This is a *wheel momentum / total momentum error* style model.
%   - All vectors are expressed in one inertial/ECI frame (simplified: body=inertial).

if nargin < 2 || strlength(ctrlMode) == 0
    ctrlMode = "baseline";
end
ctrlMode = string(ctrlMode);

% Pre-allocate
x = zeros(3, P.N);     % momentum error state
m = zeros(3, P.N);     % dipole command
tau_ext_log = zeros(3, P.N);
tau_mtq_log = zeros(3, P.N);
tau_total_log = zeros(3, P.N);
B_log = zeros(3, P.N);

solve_time_s = zeros(1, P.N);

% Initial condition
if isfield(P, "x0") && numel(P.x0) == 3
    x(:,1) = P.x0(:);
else
    x(:,1) = [0.02; -0.01; 0.03];
end

% Providers
B_provider   = @(kk) orbit_propagator_simple(kk, P);
tau_provider = @(kk) ext_torques(kk, P);

m_prev = zeros(3,1);
U_prev = [];
m_hold = zeros(3,1);

update_every = 1;
warm_start = false;
if isfield(P,'mpc') && isfield(P.mpc,'update_every')
    update_every = max(1, round(P.mpc.update_every));
end
if isfield(P,'mpc') && isfield(P.mpc,'warm_start')
    warm_start = logical(P.mpc.warm_start);
end

for k = 1:P.N
    Bk = B_provider(k);
    tk = tau_provider(k);

    B_log(:,k) = Bk;
    tau_ext_log(:,k) = tk;

    if ctrlMode == "baseline"
        mk = baseline_mtq(x(:,k), Bk, P);
    elseif ctrlMode == "mpc"
        doSolve = (mod(k-1, update_every) == 0) || (k == 1);
        if doSolve
            t0 = tic;
            U0 = [];
            if warm_start && ~isempty(U_prev)
                % Shift previous solution forward by one input and repeat the last input.
                U0 = [U_prev(4:end); U_prev(end-2:end)];
            end
            [mk, U_opt] = mpc_mtq_qp(x(:,k), k, P, B_provider, tau_provider, m_prev, U0);
            solve_time_s(k) = toc(t0);
            U_prev = U_opt;
            m_hold = mk;
        else
            mk = m_hold;
        end
    else
        error("Unknown ctrlMode '%s'. Use 'baseline' or 'mpc'.", ctrlMode);
    end

    m(:,k) = mk;
    tau_mtq_log(:,k) = cross(mk, Bk);
    tau_total_log(:,k) = tk + tau_mtq_log(:,k);

    if k < P.N
        x(:,k+1) = x(:,k) + P.Ts*( tk - skew(Bk)*mk );
    end

    m_prev = mk;
end

t = (0:P.N-1) * P.Ts;

S = struct();
S.P = P;
S.ctrlMode = ctrlMode;
S.t = t;
S.x = x;
S.m = m;
S.B_log = B_log;
S.tau_ext_log = tau_ext_log;
S.tau_mtq_log = tau_mtq_log;
S.tau_total_log = tau_total_log;
S.solve_time_s = solve_time_s;
end
