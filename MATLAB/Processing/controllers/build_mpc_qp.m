function [H, f, lb, ub] = build_mpc_qp(x0, Bseq, tauseq, P, m_prev)
%BUILD_MPC_QP Build condensed QP for LTV system:
% x_{i+1} = x_i + Ts*(tau_ext_i - [B_i]_x * m_i)
% with cost sum x'Qx + m'Rm + (dm)'S(dm)

Nh = P.mpc.Nh;
Ts = P.Ts;

Q = P.mpc.Q; R = P.mpc.R; S = P.mpc.S;
useDelta = P.mpc.useDelta;

nx = 3; nu = 3;

% Precompute G_i and d_i
G = zeros(nx, nu, Nh);
d = zeros(nx, Nh);

for i = 1:Nh
    Bi = Bseq(:,i);
    G(:,:,i) = -Ts * skew(Bi);      % torque mapping
    d(:,i)   =  Ts * tauseq(:,i);
end

% Build prediction: X = Abar*x0 + Bbar*U + Dbar
% Here A = I, so Abar is stacking identities
Abar = zeros(nx*Nh, nx);
Bbar = zeros(nx*Nh, nu*Nh);
Dbar = zeros(nx*Nh, 1);

xk = x0(:);

for i = 1:Nh
    % x_i = x0 + sum_{j=1..i} (G_j m_j + d_j)
    Abar((i-1)*nx+(1:nx), :) = eye(nx);

    % Bbar row block:
    for j = 1:i
        Bbar((i-1)*nx+(1:nx), (j-1)*nu+(1:nu)) = G(:,:,j);
    end

    % Dbar:
    Dbar((i-1)*nx+(1:nx), :) = sum(d(:,1:i), 2);
end

% Cost on states: X'Qbar X
Qbar = kron(eye(Nh), Q);
Rbar = kron(eye(Nh), R);

% Delta-m penalty
if useDelta
    % Dm = T*U - [m_prev; 0; ...]
    Tdm = zeros(nu*(Nh-1), nu*Nh);
    for i = 1:Nh-1
        Tdm((i-1)*nu+(1:nu), (i-1)*nu+(1:nu)) = -eye(nu);
        Tdm((i-1)*nu+(1:nu), (i)*nu+(1:nu))   =  eye(nu);
    end
    Sbar = kron(eye(Nh-1), S);

    % (m1 - m_prev) term
    dm0 = zeros(nu*(Nh-1),1);
    dm0(1:nu) = -m_prev(:);
else
    Tdm = [];
    Sbar = [];
    dm0 = [];
end

% Condense to QP: 1/2 U' H U + f' U
H = 2*(Bbar' * Qbar * Bbar + Rbar);
f = 2*(Bbar' * Qbar * (Abar*xk + Dbar));

if useDelta
    H = H + 2*(Tdm' * Sbar * Tdm);
    f = f + 2*(-Tdm' * Sbar * dm0);
end

% Bounds
lb = repmat(-P.m_max(:), Nh, 1);
ub = repmat( P.m_max(:), Nh, 1);

% Make symmetric for quadprog
H = 0.5*(H + H');
end
