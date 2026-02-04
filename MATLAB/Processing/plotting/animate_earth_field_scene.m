function fig = animate_earth_field_scene(P, t_s, x, m, B_T, label)
%ANIMATE_EARTH_FIELD_SCENE Interactive 3D animation: Earth + B-field + satellite + vectors.
%
% Features:
%   - 36 orbit markers (configurable via P.viz.nOrbitMarkers)
%   - B-field direction arrows along the orbit
%   - Moving "satellite" marker with live magnitude readout
%   - Play/Pause, Stop, and frame slider
%
% All vectors are shown in a single inertial/ECI frame (body = inertial in this simplified model).

if nargin < 6
    label = "run";
end

% Normalize inputs
x = ensure_nx3(x);
m = ensure_nx3(m);
B_T = ensure_nx3(B_T);
t_s = t_s(:);

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

% Keep marker/arrow count small for interactive performance.
idx_mark = unique(round(linspace(1, N, nMarkers)));

% Precompute orbit in ECI
r_eci = zeros(N,3);
for k = 1:N
    r_eci(k,:) = orbit_state_simple(k, P).';
end

Re = P.earth.Re_m;

% --- Figure / axes ---
fig = figure('Color','w', 'Name', sprintf('MTQ 3D Animation (%s)', label));
fig.Position(3:4) = [1400 900];
ax = axes(fig); %#ok<LAXES>
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
view(ax, 35, 20);
axis(ax,'vis3d');
daspect(ax,[1 1 1]);
camproj(ax,'perspective');
if isfield(P,'viz') && isfield(P.viz,'interactive') && P.viz.interactive
    try, set(fig, 'Renderer', 'opengl'); end
    try, rotate3d(fig,'on'); end
    if isfield(P.viz,'zoom_speed')
        zspeed = double(P.viz.zoom_speed);
    else
        zspeed = 0.12;
    end
    fig.WindowScrollWheelFcn = @(~,evt) camzoom(ax, (1+zspeed) ^ (-evt.VerticalScrollCount));
    try
        ax.Toolbar.Visible = 'on';
        axtoolbar(ax, {'rotate','pan','zoomin','zoomout','restoreview'});
    catch
        try, cameratoolbar(fig,'Show'); end
    end
end

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
        'FaceColor',[0.95 0.7 0.1], 'FaceAlpha',0.9, 'EdgeColor','none', ...
        'FaceLighting','gouraud', 'SpecularStrength',0.3, 'HandleVisibility','off');
    quiver3(ax, 0,0,0, sun_pos(1), sun_pos(2), sun_pos(3), 0, ...
        'Color',[0.95 0.7 0.1], 'LineWidth',1.1, 'MaxHeadSize',0.6, ...
        'HandleVisibility','off');
end

% Earth sphere (textured if available)
earthSurf = draw_earth_sphere(ax, Re, P);

% Dipole field lines + direction arrows
plot_dipole_field_lines(ax, Re, P);

% Orbit line + markers
plot3(ax, r_eci(:,1), r_eci(:,2), r_eci(:,3), 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
plot3(ax, r_eci(idx_mark,1), r_eci(idx_mark,2), r_eci(idx_mark,3), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor',[1 0.9 0.2], 'HandleVisibility','off');

% B-field direction arrows along the orbit markers
Lb = 0.12 * Re;
for ii = idx_mark(:).'
    [Bv, okB] = unit_vec(B_T(ii,:));
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

% Moving satellite marker
sat = scatter3(ax, r_eci(idx_anim(1),1), r_eci(idx_anim(1),2), r_eci(idx_anim(1),3), 160, 's', ...
    'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerEdgeColor','w', 'LineWidth', 0.9, ...
    'DisplayName','satellite');

% Vector handles (only these appear in the legend)
Lvec = 0.25 * Re;
[qB, qx, qm, qtau, qtauPerp, qtauPar] = init_vector_quivers(ax);
legend(ax, 'Location','northeastoutside');

% Lighting (optional)
if isfield(P,'viz') && isfield(P.viz,'use_lighting') && P.viz.use_lighting
    material(earthSurf, 'dull');
    try
        if exist('sun_pos','var')
            camlight(ax, sun_pos(1), sun_pos(2), sun_pos(3));
        else
            camlight(ax, 'headlight');
        end
        lighting(ax, 'gouraud');
    catch
    end
end

% HUD annotation with magnitudes
hud = annotation(fig,'textbox',[0.02 0.02 0.52 0.20], 'String', '', ...
    'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], ...
    'Interpreter','none');

sgtitle(sprintf('Earth Magnetic Field + MTQ Dipole / Torque Animation (%s)', label), 'FontWeight','bold');
xlabel(ax,'X [m]'); ylabel(ax,'Y [m]'); zlabel(ax,'Z [m]');

% --- Controls ---
frameIdx = 1;
speedFactor = 1; % frames per tick

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
    'Position',[1150 20 210 28], 'BackgroundColor','w', 'HorizontalAlignment','left');

% Timer drives the animation when "Play" is pressed.
tmr = timer('ExecutionMode','fixedRate', 'Period', 0.08, 'TimerFcn', @onTick);
start(tmr);

fig.CloseRequestFcn = @onClose;

% Initial draw
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
    if ~ishandle(fig)
        return;
    end
    if btnPlay.Value ~= 1
        return;
    end

    step = max(1, round(speedFactor));
    frameIdx = frameIdx + step;
    if frameIdx > nFrames
        frameIdx = 1; % loop
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
    r = r_eci(k,:);

    sat.XData = r(1);
    sat.YData = r(2);
    sat.ZData = r(3);

    % Controller vectors
    xk = x(k,:);
    mk = m(k,:);
    Bk = B_T(k,:);
    tau_mtq = cross(mk, Bk);
    tau_des = (-P.KH * xk(:)).';

    [bhat, okb] = unit_vec(Bk);
    if okb
        tau_par = dot(tau_des, bhat.') * bhat.'; % 1x3
    else
        tau_par = [0 0 0];
    end
    tau_perp = tau_des - tau_par;

    update_quiver(qB, r, Bk, Lvec, [0.0 0.45 0.74], 'B');
    update_quiver(qx, r, xk, Lvec, [0.35 0.35 0.35], 'x (momentum error)');
    update_quiver(qm, r, mk, Lvec, [0.85 0.33 0.10], 'm (dipole)');
    update_quiver(qtau, r, tau_mtq, Lvec, [0.47 0.67 0.19], '\tau_{mtq} = m×B');
    update_quiver(qtauPerp, r, tau_perp, Lvec, [0 0 0], '\tau_{des,\perp B}', '--');
    update_quiver(qtauPar, r, tau_par, Lvec, [0.7 0 0.7], '\tau_{des,\parallel B}', ':');

    % Live magnitudes
    xmag = norm(xk);
    Bmag_uT = 1e6 * norm(Bk);
    mmag = norm(mk);
    taumag_uNm = 1e6 * norm(tau_mtq);

    hud.String = sprintf([ ...
        'Frame %d/%d  (k=%d)\n' ...
        't = %.1f s\n' ...
        '||x|| = %.3e  [N·m·s]\n' ...
        'B = [%.1f %.1f %.1f] µT,  ||B|| = %.1f µT\n' ...
        '||m|| = %.3f  [A·m^2]\n' ...
        'tau_{mtq} = [%.2f %.2f %.2f] µN·m,  ||tau_{mtq}|| = %.2f µN·m\n' ...
        '\nNote: All vectors are in ECI/inertial (simplified).'], ...
        iFrame, nFrames, k, t_s(k), xmag, ...
        1e6*Bk(1), 1e6*Bk(2), 1e6*Bk(3), Bmag_uT, ...
        mmag, 1e6*tau_mtq(1), 1e6*tau_mtq(2), 1e6*tau_mtq(3), taumag_uNm);

    txtTime.String = sprintf('t = %.1f s', t_s(k));

    drawnow limitrate;
end
end

function [qB, qx, qm, qtau, qtauPerp, qtauPar] = init_vector_quivers(ax)
% Initialize quiver objects with legend entries.

qB = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.0 0.45 0.74], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','B');
qx = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.35 0.35 0.35], 'LineWidth', 1.2, 'MaxHeadSize', 0.7, 'DisplayName','x (momentum error)');
qm = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.85 0.33 0.10], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','m (dipole)');
qtau = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.47 0.67 0.19], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, 'DisplayName','\tau_{mtq} = m×B');
qtauPerp = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0 0 0], 'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.7, 'DisplayName','\tau_{des,\perp B}');
qtauPar = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color',[0.7 0 0.7], 'LineWidth', 1.2, 'LineStyle',':', 'MaxHeadSize', 0.7, 'DisplayName','\tau_{des,\parallel B}');

qB.AutoScale = 'off';
qx.AutoScale = 'off';
qm.AutoScale = 'off';
qtau.AutoScale = 'off';
qtauPerp.AutoScale = 'off';
qtauPar.AutoScale = 'off';
end

function update_quiver(q, r, v, L, color, name, lineStyle)
% Update quiver object to originate at r and point in direction of v (normalized).

if nargin < 7
    lineStyle = '-';
end

[uv, ok] = unit_vec(v);
if ~ok
    uv = [0;0;0];
end

q.XData = r(1);
q.YData = r(2);
q.ZData = r(3);
q.UData = L*uv(1);
q.VData = L*uv(2);
q.WData = L*uv(3);
q.Color = color;
q.DisplayName = name;
q.LineStyle = lineStyle;
end

function plot_dipole_field_lines(ax, Re, P)
% Analytic dipole field lines: r = L sin^2(theta) in Earth radii.

Ls = [1.6 2.0 3.0 4.5];
phis = deg2rad([0 60 120 180 240 300]);

theta = linspace(0.15, pi-0.15, 400); % colatitude

tilt = 0;
if isfield(P,'earth') && isfield(P.earth,'dipole_tilt_deg')
    tilt = deg2rad(P.earth.dipole_tilt_deg);
end
Ry = [cos(tilt) 0 sin(tilt); 0 1 0; -sin(tilt) 0 cos(tilt)]; % tilt about +Y (X-Z plane)

for L = Ls
    r = L * (sin(theta).^2);
    for phi = phis
        x = r .* sin(theta) * cos(phi);
        y = r .* sin(theta) * sin(phi);
        z = r .* cos(theta);
        pts = Ry * [Re*x; Re*y; Re*z];
        plot3(ax, pts(1,:), pts(2,:), pts(3,:), '-', 'Color',[0.2 0.7 0.9], 'LineWidth', 0.6, 'HandleVisibility','off');

        % Direction arrowheads along the line
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

% Axis indicator
quiver3(ax, 0,0,0, 1.2*Re,0,0, 0, 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.4, 'HandleVisibility','off');
quiver3(ax, 0,0,0, 0,1.2*Re,0, 0, 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.4, 'HandleVisibility','off');
quiver3(ax, 0,0,0, 0,0,1.2*Re, 0, 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.4, 'HandleVisibility','off');
text(ax, 1.25*Re,0,0,'X', 'HandleVisibility','off');
text(ax, 0,1.25*Re,0,'Y', 'HandleVisibility','off');
text(ax, 0,0,1.25*Re,'Z', 'HandleVisibility','off');

% Dipole axis indicator (tilted)
m_hat = [sin(tilt); 0; cos(tilt)];
quiver3(ax, 0,0,0, 1.2*Re*m_hat(1), 1.2*Re*m_hat(2), 1.2*Re*m_hat(3), 0, ...
    'Color',[0.2 0.2 0.2], 'LineWidth', 1.0, 'MaxHeadSize', 0.4, 'HandleVisibility','off');
text(ax, 1.25*Re*m_hat(1), 1.25*Re*m_hat(2), 1.25*Re*m_hat(3), 'dipole axis', 'HandleVisibility','off');
end

function [u, ok] = unit_vec(v)
% Normalize vector; returns ok=false for near-zero.
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
