function m_cmd = mpc_mtq_qp(xk, k, P, B_provider, tau_provider, m_prev)
%MPC_MTQ_QP Compute MPC command at step k using predicted B and tau sequences

Nh = P.mpc.Nh;

Bseq   = zeros(3, Nh);
tauseq = zeros(3, Nh);

for i = 1:Nh
    kk = k + (i-1);
    Bseq(:,i)   = B_provider(kk);
    tauseq(:,i) = tau_provider(kk);
end

[H, f, lb, ub] = build_mpc_qp(xk, Bseq, tauseq, P, m_prev);

% Solve QP
U = quadprog(H, f, [], [], [], [], lb, ub, [], P.mpc.qp_opts);

% Fallback if solver fails
if isempty(U)
    m_cmd = zeros(3,1);
else
    m_cmd = U(1:3);
end

end
