function [m_cmd, U] = mpc_mtq_qp(xk, k, P, B_provider, tau_provider, m_prev, U0)
%MPC_MTQ_QP Compute MPC command at step k using predicted B and tau sequences.
%
% Optional:
%   U0 [3*Nh x 1] initial guess for quadprog (warm-start)

Nh = P.mpc.Nh;

Bseq   = zeros(3, Nh);
tauseq = zeros(3, Nh);

for i = 1:Nh
    kk = k + (i-1);
    Bseq(:,i)   = B_provider(kk);
    tauseq(:,i) = tau_provider(kk);
end

[H, f, lb, ub] = build_mpc_qp(xk, Bseq, tauseq, P, m_prev);

% Solve QP (optionally warm-start)
if nargin < 7
    U0 = [];
end
U = quadprog(H, f, [], [], [], [], lb, ub, U0, P.mpc.qp_opts);

% Fallback if solver fails
if isempty(U)
    m_cmd = zeros(3,1);
    U = [];
else
    m_cmd = U(1:3);
end

% Optional post-processing: adjust *only* the component along B to reduce wasted dipole
% while keeping the exact same torque tau = m×B (adding alpha*B does not change tau).
if isfield(P,'mpc') && isfield(P.mpc,'enforce_tau_perpB') && P.mpc.enforce_tau_perpB
    Bk = Bseq(:,1);
    m_cmd = shift_along_B_into_box(m_cmd, Bk, P.m_max(:));
end

end

function m_adj = shift_along_B_into_box(m, B, m_max)
%SHIFT_ALONG_B_INTO_BOX Find m_adj = m_perp + alpha*B such that:
%  - cross(m_adj, B) == cross(m, B) (same torque)
%  - -m_max <= m_adj <= m_max (component-wise)
%  - alpha is chosen as close to 0 as possible (minimizes ||m_parallel||)

m = m(:);
B = B(:);
m_max = m_max(:);

Bn2 = max(dot(B,B), 1e-12);
beta = dot(m, B) / Bn2;
m_perp = m - beta * B;

% Solve bounds for alpha in each component:
%   -m_max_i <= m_perp_i + alpha*B_i <= m_max_i
amin = -inf;
amax = inf;
epsB = 1e-15;
for i = 1:3
    Bi = B(i);
    if abs(Bi) < epsB
        % This component is not affected by alpha; must already be within bounds.
        if m_perp(i) < -m_max(i) - 1e-12 || m_perp(i) > m_max(i) + 1e-12
            % Should never happen if original m was feasible; fallback.
            m_adj = min(max(m, -m_max), m_max);
            return;
        end
        continue;
    end

    a1 = (-m_max(i) - m_perp(i)) / Bi;
    a2 = ( m_max(i) - m_perp(i)) / Bi;
    lo = min(a1, a2);
    hi = max(a1, a2);
    amin = max(amin, lo);
    amax = min(amax, hi);
end

% The original alpha=beta always produces the original feasible m, so intersection should exist.
if ~(amin <= amax)
    m_adj = min(max(m, -m_max), m_max);
    return;
end

alpha = min(max(0, amin), amax); % closest to 0 in [amin,amax]
m_adj = m_perp + alpha * B;

% Final safety clamp for tiny numeric spill; should be no-op.
m_adj = min(max(m_adj, -m_max), m_max);
end
