function fig = plot_earth_field_scene(P, t_s, x, m, B_T, label)
%PLOT_EARTH_FIELD_SCENE 3D scene: Earth, dipole field lines, orbit, and MTQ vectors.
%
% This is a visualization aid. Vectors are directionally correct; arrow lengths
% are normalized for readability.

if nargin < 6
    label = "run";
end

% Normalize
x = ensure_nx3(x);
m = ensure_nx3(m);
B_T = ensure_nx3(B_T);
t_s = t_s(:);

N = numel(t_s);

% Use a representative snapshot (start and end)
k0 = 1;
kf = N;

% Orbit marker count (user-requested default: 36)
nMarkers = 36;
if isfield(P, "viz") && isfield(P.viz, "nOrbitMarkers")
    nMarkers = max(2, round(P.viz.nOrbitMarkers));
end

% Orbit curve
r_eci = zeros(N,3);
for k = 1:N
    r_eci(k,:) = orbit_state_simple(k, P).';
end

Re = P.earth.Re_m;

fig = figure('Color','w','Name',sprintf('Earth / B-field / MTQ Scene (%s)', label));
fig.Position(3:4) = [1400 900];
ax = axes(fig); %#ok<LAXES>
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
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

% Click-to-focus: set camera target to the point under the cursor.
try
    fig.WindowButtonDownFcn = @(~,~) set_camtarget_to_cursor(ax);
catch
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

% A few ideal dipole field lines (analytic r = L sin^2(theta))
plot_dipole_field_lines(ax, Re, P);

% Orbit
plot3(ax, r_eci(:,1), r_eci(:,2), r_eci(:,3), 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
idx_mark = unique(round(linspace(1, N, nMarkers)));
nMarkersUsed = numel(idx_mark);
plot3(ax, r_eci(idx_mark,1), r_eci(idx_mark,2), r_eci(idx_mark,3), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor',[1 0.9 0.2], 'HandleVisibility','off');
scatter3(ax, r_eci(k0,1), r_eci(k0,2), r_eci(k0,3), 90, 's', 'MarkerFaceColor','g', ...
    'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');
scatter3(ax, r_eci(kf,1), r_eci(kf,2), r_eci(kf,3), 90, 's', 'MarkerFaceColor','r', ...
    'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');

% B-field direction arrows along the orbit (flow direction at the satellite)
Lb = 0.12 * Re;
for ii = idx_mark(:).'
    [Bv, okB] = unit_vec(B_T(ii,:));
    if okB
        quiver3(ax, r_eci(ii,1), r_eci(ii,2), r_eci(ii,3), Lb*Bv(1), Lb*Bv(2), Lb*Bv(3), 0, ...
            'Color',[0.0 0.45 0.74], 'LineWidth', 0.7, 'MaxHeadSize', 0.6, 'HandleVisibility','off');
    end
end

% Vectors at start and end (only the start marker contributes legend entries)
plot_vectors(ax, r_eci(k0,:), x(k0,:), m(k0,:), B_T(k0,:), P, 'start');
plot_vectors(ax, r_eci(kf,:), x(kf,:), m(kf,:), B_T(kf,:), P, 'end');

xlabel(ax,'X [m]'); ylabel(ax,'Y [m]'); zlabel(ax,'Z [m]');
view(ax, 35, 20);

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

x0n = norm(x(k0,:)); xfn = norm(x(kf,:));

sgtitle(sprintf('Earth Magnetic Field + MTQ Dipole / Torque Vectors (%s)', label), 'FontWeight','bold');

% Text box with key magnitudes
B0_uT = 1e6*norm(B_T(k0,:));
Bf_uT = 1e6*norm(B_T(kf,:));

m0 = norm(m(k0,:));
mf = norm(m(kf,:));

tau0_uNm = 1e6*norm(cross(m(k0,:), B_T(k0,:)));
tauf_uNm = 1e6*norm(cross(m(kf,:), B_T(kf,:)));

txt = sprintf([ ...
    'Markers: start (green) and end (red)\n' ...
    'Orbit markers: %d evenly spaced samples\n' ...
    'x magnitude: start %.3e, end %.3e [N·m·s]\n' ...
    '|B|: start %.1f, end %.1f [µT]\n' ...
    '|m|: start %.3f, end %.3f [A·m^2]\n' ...
    '|tau_{mtq}|: start %.2f, end %.2f [µN·m]\n' ...
    '\nVectors shown at each marker:\n' ...
    '  B (blue), x (gray), m (red), tau_{mtq}=m×B (green)\n' ...
    '  tau_{des,⊥B} (black dashed) and tau_{des,∥B} (magenta dotted)\n' ...
    '\nBaseline note: MTQ torque is always ⟂B, so the ∥B component is not directly achievable.\n' ...
    'All vectors are shown in a single inertial/ECI frame (body = inertial in this simplified model).\n' ...
    'Arrow lengths are normalized for readability.' ...
    ], nMarkersUsed, x0n, xfn, B0_uT, Bf_uT, m0, mf, tau0_uNm, tauf_uNm);
annotation(fig,'textbox',[0.02 0.02 0.55 0.22], 'String', txt, ...
    'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], ...
    'Interpreter','none');
end

function plot_vectors(ax, r0, x0, m0, B0, P, tag)
% Plot normalized direction vectors at point r0.

KH = P.KH;
tau_des = (-KH * x0(:)).';
tau_mtq = cross(m0, B0);

[bhat, okb] = unit_vec(B0);
if okb
    tau_par = dot(tau_des, bhat.') * bhat.'; % 1x3
else
    tau_par = [0 0 0];
end
tau_perp = tau_des - tau_par;

% Normalize arrow lengths for readability
Re = P.earth.Re_m;
L = 0.25 * Re;

[Bv, okB] = unit_vec(B0);
[mv, okm] = unit_vec(m0);
[tv, okt] = unit_vec(tau_mtq);
[dpv, okdp] = unit_vec(tau_perp);
[dav, okda] = unit_vec(tau_par);
[xv, okx] = unit_vec(x0);

showLegend = strcmp(tag,'start');
hv = 'off';
if showLegend
    hv = 'on';
end

if okB
    h = quiver3(ax, r0(1), r0(2), r0(3), L*Bv(1), L*Bv(2), L*Bv(3), 0, ...
        'Color',[0.0 0.45 0.74], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','B');
    try, h.Clipping = 'off'; end
end
if okx
    h = quiver3(ax, r0(1), r0(2), r0(3), L*xv(1), L*xv(2), L*xv(3), 0, ...
        'Color',[0.35 0.35 0.35], 'LineWidth', 1.2, 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','x (momentum error)');
    try, h.Clipping = 'off'; end
end
if okm
    h = quiver3(ax, r0(1), r0(2), r0(3), L*mv(1), L*mv(2), L*mv(3), 0, ...
        'Color',[0.85 0.33 0.10], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','m (dipole)');
    try, h.Clipping = 'off'; end
end
if okt
    h = quiver3(ax, r0(1), r0(2), r0(3), L*tv(1), L*tv(2), L*tv(3), 0, ...
        'Color',[0.47 0.67 0.19], 'LineWidth', 1.3, 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','\tau_{mtq} = m×B');
    try, h.Clipping = 'off'; end
end
if okdp
    h = quiver3(ax, r0(1), r0(2), r0(3), L*dpv(1), L*dpv(2), L*dpv(3), 0, ...
        'Color',[0 0 0], 'LineWidth', 1.2, 'LineStyle','--', 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','\tau_{des,\perp B}');
    try, h.Clipping = 'off'; end
end
if okda
    h = quiver3(ax, r0(1), r0(2), r0(3), L*dav(1), L*dav(2), L*dav(3), 0, ...
        'Color',[0.7 0 0.7], 'LineWidth', 1.2, 'LineStyle',':', 'MaxHeadSize', 0.7, ...
        'HandleVisibility', hv, 'DisplayName','\tau_{des,\parallel B}');
    try, h.Clipping = 'off'; end
end

text(ax, r0(1), r0(2), r0(3), ['  ' tag], 'Color','k', 'FontWeight','bold');
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

        % Add a few arrowheads along the line to indicate field "flow" direction.
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

function set_camtarget_to_cursor(ax)
try
    cp = ax.CurrentPoint;
    ax.CameraTarget = cp(1,1:3);
catch
end
end
