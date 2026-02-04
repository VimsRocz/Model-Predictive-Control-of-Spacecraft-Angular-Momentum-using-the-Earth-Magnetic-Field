function fig = animate_mtq_full_scene(P, S, label)
%ANIMATE_MTQ_FULL_SCENE Interactive 3D animation for the full model.
%
% Shows:
% - Earth + dipole field lines (context) + optional Sun
% - Orbit markers or full-rotation animation (configurable)
% - Satellite cube + body axes
% - Vectors: B, h_w, m, tau_mtq, tau_des, and the plane ⟂B
% - Magnetometer component arrows (B_body on body axes)
% - Live HUD: magnetometer components, torque components + totals, tau·B (≈0)

if nargin < 3
    label = "full";
end

t_s = S.t(:);
r_eci = ensure_nx3(S.r_eci);
B_eci = ensure_nx3(S.B_eci);

N = numel(t_s);

% Orbit sample selection (full rotation or markers)
nMarkers = 36;
if isfield(P, "viz") && isfield(P.viz, "nOrbitMarkers")
    nMarkers = max(2, round(P.viz.nOrbitMarkers));
end
use_all = isfield(P,"viz") && isfield(P.viz,"animate_all_steps") && P.viz.animate_all_steps;
stride = 1;
if isfield(P,"viz") && isfield(P.viz,"animate_stride")
    stride = max(1, round(P.viz.animate_stride));
end
if use_all
    idx_anim = 1:stride:N;
else
    idx_anim = unique(round(linspace(1, N, nMarkers)));
end
nFrames = numel(idx_anim);

Re = P.earth.Re_m;

% --- Figure / axes ---
fig = figure('Color','w', 'Name', sprintf('Full MTQ 3D Animation (%s)', label));
fig.Position(3:4) = [1500 980];
ax = axes(fig); %#ok<LAXES>
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
view(ax, 35, 20);

% Sun (optional)
if isfield(P,'viz') && isfield(P.viz,'show_sun') && P.viz.show_sun
    sun_dir = [1;0;0];
    if isfield(P.viz,'sun_dir_eci') && numel(P.viz.sun_dir_eci) == 3
        sun_dir = P.viz.sun_dir_eci(:);
    end
    sun_dir = sun_dir / max(norm(sun_dir), 1e-12);
    sun_dist = 20;
    if isfield(P.viz,'sun_distance_Re'); sun_dist = double(P.viz.sun_distance_Re); end
    sun_rad = 2.5;
    if isfield(P.viz,'sun_radius_Re'); sun_rad = double(P.viz.sun_radius_Re); end
    sun_pos = (sun_dist * Re) * sun_dir;
    [xs, ys, zs] = sphere(40);
    surf(ax, sun_pos(1)+sun_rad*Re*xs, sun_pos(2)+sun_rad*Re*ys, sun_pos(3)+sun_rad*Re*zs, ...
        'FaceColor',[0.95 0.7 0.1], 'FaceAlpha',0.85, 'EdgeColor','none', ...
        'HandleVisibility','off');
    quiver3(ax, 0,0,0, sun_pos(1), sun_pos(2), sun_pos(3), 0, ...
        'Color',[0.95 0.7 0.1], 'LineWidth',1.1, 'MaxHeadSize',0.6, ...
        'HandleVisibility','off');
end

% Earth sphere
[xe, ye, ze] = sphere(60);
surf(ax, Re*xe, Re*ye, Re*ze, ...
    'FaceColor',[0.15 0.35 0.7], 'FaceAlpha',0.25, 'EdgeColor','none', ...
    'HandleVisibility','off');

% Dipole field lines + arrows (context)
plot_dipole_field_lines(ax, Re, P);

% Orbit line + markers
plot3(ax, r_eci(:,1), r_eci(:,2), r_eci(:,3), 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
plot3(ax, r_eci(idx_anim,1), r_eci(idx_anim,2), r_eci(idx_anim,3), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor',[1 0.9 0.2], 'HandleVisibility','off');

% B direction arrows along orbit markers
Lb = 0.12 * Re;
for ii = idx_anim(:).'
    [Bv, okB] = unit_vec(B_eci(ii,:));
    if okB
        quiver3(ax, r_eci(ii,1), r_eci(ii,2), r_eci(ii,3), Lb*Bv(1), Lb*Bv(2), Lb*Bv(3), 0, ...
            'Color',[0.0 0.45 0.74], 'LineWidth', 0.7, 'MaxHeadSize', 0.6, 'HandleVisibility','off');
    end
end

% Start/end markers
scatter3(ax, r_eci(idx_anim(1),1), r_eci(idx_anim(1),2), r_eci(idx_anim(1),3), 90, 's', ...
    'MarkerFaceColor','g', 'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');
scatter3(ax, r_eci(idx_anim(end),1), r_eci(idx_anim(end),2), r_eci(idx_anim(end),3), 90, 's', ...
    'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');

% --- Satellite cube (patch) ---
satSizeRe = 0.02;
if isfield(P,'viz') && isfield(P.viz,'sat_size_Re')
    satSizeRe = double(P.viz.sat_size_Re);
end
satSize_m = Re * satSizeRe;

[Vb, F] = unit_cube();
satPatch = patch(ax, 'Vertices', zeros(8,3), 'Faces', F, ...
    'FaceColor',[0.2 0.2 0.2], 'FaceAlpha', 0.85, 'EdgeColor',[1 1 1]*0.9, ...
    'LineWidth', 0.4, 'DisplayName','satellite');

% Plane ⟂B patch
planePatch = patch(ax, 'XData',nan(1,4), 'YData',nan(1,4), 'ZData',nan(1,4), [0.6 0.6 0.6], ...
    'FaceAlpha', 0.12, 'EdgeColor',[0.5 0.5 0.5], 'LineStyle','--', ...
    'DisplayName','plane ⟂B');

% Body axes quivers
qbx = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.85 0.33 0.10], 'LineWidth', 1.0, 'MaxHeadSize', 0.5, 'DisplayName','body x');
qby = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.25 0.6 0.2],  'LineWidth', 1.0, 'MaxHeadSize', 0.5, 'DisplayName','body y');
qbz = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.1 0.3 0.8],   'LineWidth', 1.0, 'MaxHeadSize', 0.5, 'DisplayName','body z');
qbx.AutoScale = 'off'; qby.AutoScale = 'off'; qbz.AutoScale = 'off';

% Magnetometer component arrows (B_body along body axes)
qmx = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.85 0.33 0.10], 'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.6, 'DisplayName','B_x (mag)');
qmy = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.25 0.6 0.2],  'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.6, 'DisplayName','B_y (mag)');
qmz = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.1 0.3 0.8],   'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.6, 'DisplayName','B_z (mag)');
qmx.AutoScale = 'off'; qmy.AutoScale = 'off'; qmz.AutoScale = 'off';

% Controller vectors quivers (legend entries)
Lvec = 0.28 * Re;
qB   = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.0 0.45 0.74], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','B');
qh   = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.35 0.35 0.35], 'LineWidth', 1.2, 'MaxHeadSize', 0.7, 'DisplayName','h_w');
qm   = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.85 0.33 0.10], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','m');
qtau = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.47 0.67 0.19], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','tau_{mtq}');
qtauPerp = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0 0 0], 'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.7, 'DisplayName','tau_{des,⊥B}');
qtauPar  = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.7 0 0.7], 'LineWidth', 1.2, 'LineStyle',':',  'MaxHeadSize', 0.7, 'DisplayName','tau_{des,∥B}');
qB.AutoScale = 'off'; qh.AutoScale = 'off'; qm.AutoScale = 'off'; qtau.AutoScale = 'off'; qtauPerp.AutoScale = 'off'; qtauPar.AutoScale = 'off';

legend(ax, 'Location','northeastoutside');

sgtitle(sprintf('Full MTQ Desaturation Animation (%s)', label), 'FontWeight','bold');
xlabel(ax,'X [m]'); ylabel(ax,'Y [m]'); zlabel(ax,'Z [m]');

% HUD
hud = annotation(fig,'textbox',[0.02 0.02 0.56 0.22], 'String', '', ...
    'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], ...
    'Interpreter','none');

% --- Controls ---
frameIdx = 1;
speedFactor = 1;

btnPlay = uicontrol(fig, 'Style','togglebutton', 'String','Play', ...
    'Units','pixels', 'Position',[20 20 70 28], 'Callback', @onPlayToggle);
btnStop = uicontrol(fig, 'Style','pushbutton', 'String','Stop', ...
    'Units','pixels', 'Position',[100 20 70 28], 'Callback', @onStop);

uicontrol(fig, 'Style','text', 'String','Speed', 'Units','pixels', ...
    'Position',[180 23 45 20], 'BackgroundColor','w', 'HorizontalAlignment','left');
popupSpeed = uicontrol(fig, 'Style','popupmenu', 'String',{'0.5x','1x','2x','4x'}, ...
    'Units','pixels', 'Position',[230 20 70 28], 'Callback', @onSpeedChange, 'Value', 2);

slider = uicontrol(fig, 'Style','slider', 'Min',1, 'Max',max(2,nFrames), 'Value',1, ...
    'Units','pixels', 'Position',[320 20 820 28], 'Callback', @onSliderMove);
if nFrames > 1
    slider.SliderStep = [1/(nFrames-1) min(10, nFrames-1)/(nFrames-1)];
end
txtTime = uicontrol(fig, 'Style','text', 'String','t = 0 s', 'Units','pixels', ...
    'Position',[1150 20 320 28], 'BackgroundColor','w', 'HorizontalAlignment','left');

tmr = timer('ExecutionMode','fixedRate', 'Period', 0.08, 'TimerFcn', @onTick);
start(tmr);
fig.CloseRequestFcn = @onClose;

update_frame(1);

% ---------------- nested callbacks ----------------
function onPlayToggle(~, ~)
    if btnPlay.Value == 1
        btnPlay.String = 'Pause';
    else
        btnPlay.String = 'Play';
    end
end

function onStop(~, ~)
    btnPlay.Value = 0;
    btnPlay.String = 'Play';
    frameIdx = 1;
    slider.Value = 1;
    update_frame(frameIdx);
end

function onSpeedChange(~, ~)
    switch popupSpeed.Value
        case 1, speedFactor = 0.5;
        case 2, speedFactor = 1;
        case 3, speedFactor = 2;
        case 4, speedFactor = 4;
        otherwise, speedFactor = 1;
    end
end

function onSliderMove(~, ~)
    frameIdx = clamp_frame(round(slider.Value));
    slider.Value = frameIdx;
    update_frame(frameIdx);
end

function onTick(~, ~)
    if ~ishandle(fig); return; end
    if btnPlay.Value ~= 1; return; end

    step = max(1, round(speedFactor));
    frameIdx = frameIdx + step;
    if frameIdx > nFrames
        frameIdx = 1;
    end
    slider.Value = frameIdx;
    update_frame(frameIdx);
end

function onClose(~, ~)
    try
        stop(tmr);
        delete(tmr);
    catch
    end
    delete(fig);
end

% ---------------- helpers ----------------
function i = clamp_frame(i)
    i = max(1, min(nFrames, i));
end

function update_frame(iFrame)
    k = idx_anim(iFrame);
    r = r_eci(k,:).';
    qk = S.q_ib(k,:);
    C_bi = quat2dcm(qk);
    C_ib = C_bi.';

    % Satellite cube
    V = (C_ib * (satSize_m * Vb)) + r;
    satPatch.Vertices = V.';

    % Body axes triad
    Lbdy = 0.16 * Re;
    ex = C_ib(:,1); ey = C_ib(:,2); ez = C_ib(:,3);
    update_quiver_dir(qbx, r, ex, Lbdy);
    update_quiver_dir(qby, r, ey, Lbdy);
    update_quiver_dir(qbz, r, ez, Lbdy);

    % Vectors (convert body->ECI where needed)
    Bk = B_eci(k,:).';
    Bk_body = S.B_body(k,:).';
    mk_body = S.m(k,:).';
    mk_eci = C_ib * mk_body;
    hk_body = S.h_w(k,:).';
    hk_eci = C_ib * hk_body;
    tau_mtq_eci = C_ib * S.tau_mtq_body(k,:).';
    tau_mtq_body = S.tau_mtq_body(k,:).';
    tau_ext_body = S.tau_ext_body(k,:).';
    tau_rw_body = S.tau_rw_body(k,:).';
    tau_total_body = tau_ext_body + tau_mtq_body + tau_rw_body;

    tau_des_body = -P.KH * hk_body;
    tau_des_eci = C_ib * tau_des_body;

    bhat = Bk / max(norm(Bk), 1e-30);
    tau_par = dot(tau_des_eci, bhat) * bhat;
    tau_perp = tau_des_eci - tau_par;

    update_quiver(qB, r, Bk, Lvec);
    update_quiver(qh, r, hk_eci, Lvec);
    update_quiver(qm, r, mk_eci, Lvec);
    update_quiver(qtau, r, tau_mtq_eci, Lvec);
    update_quiver(qtauPerp, r, tau_perp, Lvec);
    update_quiver(qtauPar, r, tau_par, Lvec);

    % Plane ⟂B
    [X,Y,Z] = plane_quad(r, bhat, 0.14*Re);
    planePatch.XData = X; planePatch.YData = Y; planePatch.ZData = Z;

    % Magnetometer component arrows along body axes
    Lmag = 0.18 * Re;
    Bn = max(norm(Bk_body), 1e-30);
    bx = (Bk_body(1)/Bn) * C_ib(:,1);
    by = (Bk_body(2)/Bn) * C_ib(:,2);
    bz = (Bk_body(3)/Bn) * C_ib(:,3);
    update_quiver_dir(qmx, r, bx, Lmag);
    update_quiver_dir(qmy, r, by, Lmag);
    update_quiver_dir(qmz, r, bz, Lmag);

    % HUD values
    tau_dotB = dot(tau_mtq_eci, bhat);
    ang_deg = acosd(max(-1, min(1, dot(tau_mtq_eci, bhat) / max(norm(tau_mtq_eci), 1e-30))));
    qe = quatmultiply(quatconj(S.q_ib_des(k,:)), qk);
    if qe(1) < 0; qe = -qe; end
    att_err_deg = 2 * acosd(max(-1, min(1, qe(1))));

    hud.String = sprintf([ ...
        'Frame %d/%d  (k=%d)\n' ...
        't = %.1f s\n' ...
        '||h_w|| = %.3e [N·m·s]\n' ...
        'Magnetometer B_body = [%.1f %.1f %.1f] µT,  ||B|| = %.1f µT\n' ...
        '||m|| = %.3f [A·m^2]\n' ...
        'tau_{mtq,body} = [%.2f %.2f %.2f] µN·m,  ||tau_{mtq}|| = %.2f µN·m\n' ...
        'tau_{total,body} = [%.2f %.2f %.2f] µN·m,  ||tau_total|| = %.2f µN·m\n' ...
        'tau_{mtq}·bhat = %.2e (ideal 0)\n' ...
        'angle(tau_{mtq},B) = %.2f deg (ideal 90)\n' ...
        'attitude error = %.3f deg\n' ...
        '\nTip: use mouse/scroll to rotate + zoom the 3D scene.' ...
        ], iFrame, nFrames, k, t_s(k), norm(hk_body), ...
        1e6*Bk_body(1), 1e6*Bk_body(2), 1e6*Bk_body(3), 1e6*norm(Bk), ...
        norm(mk_body), ...
        1e6*tau_mtq_body(1), 1e6*tau_mtq_body(2), 1e6*tau_mtq_body(3), 1e6*norm(tau_mtq_body), ...
        1e6*tau_total_body(1), 1e6*tau_total_body(2), 1e6*tau_total_body(3), 1e6*norm(tau_total_body), ...
        tau_dotB, ang_deg, att_err_deg);

    txtTime.String = sprintf('t = %.1f s  |  mode=%s  |  bfield=%s', t_s(k), string(S.desatMode), bfield_label(P));

    drawnow limitrate;
end
end

function update_quiver(q, r, v, L)
[u, ok] = unit_vec(v);
if ~ok; u = [0;0;0]; end
q.XData = r(1); q.YData = r(2); q.ZData = r(3);
q.UData = L*u(1); q.VData = L*u(2); q.WData = L*u(3);
end

function update_quiver_dir(q, r, dir, L)
q.XData = r(1); q.YData = r(2); q.ZData = r(3);
q.UData = L*dir(1); q.VData = L*dir(2); q.WData = L*dir(3);
end

function [X,Y,Z] = plane_quad(r, n_hat, halfSize)
a = [0;0;1];
if abs(dot(a, n_hat)) > 0.9
    a = [0;1;0];
end
u = cross(n_hat, a);
u = u / max(norm(u), 1e-12);
v = cross(n_hat, u);
v = v / max(norm(v), 1e-12);
p1 = r + halfSize*( u + v);
p2 = r + halfSize*( u - v);
p3 = r + halfSize*(-u - v);
p4 = r + halfSize*(-u + v);
X = [p1(1) p2(1) p3(1) p4(1)];
Y = [p1(2) p2(2) p3(2) p4(2)];
Z = [p1(3) p2(3) p3(3) p4(3)];
end

function [Vb, F] = unit_cube()
% Unit cube vertices as 3x8 in body coords (centered).
Vb = [ ...
    -0.5 -0.5 -0.5;
     0.5 -0.5 -0.5;
     0.5  0.5 -0.5;
    -0.5  0.5 -0.5;
    -0.5 -0.5  0.5;
     0.5 -0.5  0.5;
     0.5  0.5  0.5;
    -0.5  0.5  0.5 ].';

F = [ ...
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8 ];
end

function s = bfield_label(P)
if isfield(P,'env') && isfield(P.env,'bfield_model')
    s = string(P.env.bfield_model);
else
    s = "dipole";
end
end

function [u, ok] = unit_vec(v)
v = v(:);
n = norm(v);
if n < 1e-30
    u = [0;0;0];
    ok = false;
else
    u = v / n;
    ok = true;
end
end

function A = ensure_nx3(A)
if size(A,2) == 3
    return;
end
if size(A,1) == 3
    A = A.';
    return;
end
error('Expected Nx3 or 3xN array. Got %dx%d.', size(A,1), size(A,2));
end

function plot_dipole_field_lines(ax, Re, P)
Ls = [1.6 2.0 3.0 4.5];
phis = deg2rad([0 60 120 180 240 300]);
theta = linspace(0.15, pi-0.15, 400);

tilt = 0;
if isfield(P,'earth') && isfield(P.earth,'dipole_tilt_deg')
    tilt = deg2rad(P.earth.dipole_tilt_deg);
end
Ry = [cos(tilt) 0 sin(tilt); 0 1 0; -sin(tilt) 0 cos(tilt)];

for L = Ls
    rr = L * (sin(theta).^2);
    for phi = phis
        x = rr .* sin(theta) * cos(phi);
        y = rr .* sin(theta) * sin(phi);
        z = rr .* cos(theta);
        pts = Ry * [Re*x; Re*y; Re*z];
        plot3(ax, pts(1,:), pts(2,:), pts(3,:), '-', 'Color',[0.2 0.7 0.9], 'LineWidth', 0.6, 'HandleVisibility','off');

        arrowLen = 0.07 * Re;
        idxA = round(linspace(80, 320, 3));
        for kk = idxA(:).'
            p = pts(:,kk);
            Bp = earth_dipole_field(p, P);
            [bv, ok] = unit_vec(Bp);
            if ok
                quiver3(ax, p(1), p(2), p(3), arrowLen*bv(1), arrowLen*bv(2), arrowLen*bv(3), 0, ...
                    'Color',[0.2 0.7 0.9], 'LineWidth', 0.5, 'MaxHeadSize', 0.55, 'HandleVisibility','off');
            end
        end
    end
end
end
