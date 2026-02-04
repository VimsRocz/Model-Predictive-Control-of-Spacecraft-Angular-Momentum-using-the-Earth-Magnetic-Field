function fig = plot_earth_field_scene_full(P, S, label)
%PLOT_EARTH_FIELD_SCENE_FULL 3D scene (full model): Earth, orbit, satellite, B, MTQ torque.
%
% Shows the "MTQ torque is always perpendicular to B" limitation using a
% plane normal to B at the satellite position, plus key vectors.

if nargin < 3
    label = "full";
end

t = S.t(:);
r_eci = ensure_nx3(S.r_eci);

N = size(r_eci,1);
k0 = 1;
kf = N;

Re = P.earth.Re_m;

% Orbit marker count
nMarkers = 36;
if isfield(P, "viz") && isfield(P.viz, "nOrbitMarkers")
    nMarkers = max(2, round(P.viz.nOrbitMarkers));
end
idx_mark = unique(round(linspace(1, N, nMarkers)));

fig = figure('Color','w','Name',sprintf('Full 3D Scene (%s)', label));
fig.Position(3:4) = [1450 950];
ax = axes(fig); %#ok<LAXES>
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
view(ax, 35, 20);
axis(ax,'vis3d');
daspect(ax,[1 1 1]);
camproj(ax,'perspective');

% Interactive helpers
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

% Dipole field lines for context (even if IGRF is used for simulation)
plot_dipole_field_lines(ax, Re, P);

% Orbit
plot3(ax, r_eci(:,1), r_eci(:,2), r_eci(:,3), 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
plot3(ax, r_eci(idx_mark,1), r_eci(idx_mark,2), r_eci(idx_mark,3), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor',[1 0.9 0.2], 'HandleVisibility','off');

scatter3(ax, r_eci(k0,1), r_eci(k0,2), r_eci(k0,3), 90, 's', 'MarkerFaceColor','g', ...
    'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');
scatter3(ax, r_eci(kf,1), r_eci(kf,2), r_eci(kf,3), 90, 's', 'MarkerFaceColor','r', ...
    'MarkerEdgeColor','k', 'LineWidth', 0.6, 'HandleVisibility','off');

% B-field direction arrows at markers (from simulation, ECI)
if isfield(S,'B_eci')
    B_eci = ensure_nx3(S.B_eci);
else
    B_eci = zeros(N,3);
end

Lb = 0.12 * Re;
for ii = idx_mark(:).'
    [Bv, okB] = unit_vec(B_eci(ii,:));
    if okB
        quiver3(ax, r_eci(ii,1), r_eci(ii,2), r_eci(ii,3), Lb*Bv(1), Lb*Bv(2), Lb*Bv(3), 0, ...
            'Color',[0.0 0.45 0.74], 'LineWidth', 0.7, 'MaxHeadSize', 0.6, 'HandleVisibility','off');
    end
end

% Vectors at start and end
Lvec = 0.28 * Re;
[~, ~, txt0] = plot_vectors_at_k(ax, P, S, k0, Lvec, true);
[~, ~, txtf] = plot_vectors_at_k(ax, P, S, kf, Lvec, false);

% Satellite cubes at start/end
satSizeRe = 0.02;
if isfield(P,'viz') && isfield(P.viz,'sat_size_Re')
    satSizeRe = double(P.viz.sat_size_Re);
end
draw_sat_cube(ax, r_eci(k0,:).', S.q_ib(k0,:), Re*satSizeRe, [0.2 0.2 0.2], 0.9, "sat (start)");
draw_sat_cube(ax, r_eci(kf,:).', S.q_ib(kf,:), Re*satSizeRe, [0.2 0.2 0.2], 0.6, "sat (end)");

xlabel(ax,'X [m]'); ylabel(ax,'Y [m]'); zlabel(ax,'Z [m]');
legend(ax, 'Location','northeastoutside');
sgtitle(sprintf('Earth Magnetic Field + MTQ Desaturation (full model) (%s)', label), 'FontWeight','bold');

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

% HUD textbox
txt = sprintf([ ...
    'Markers: start (green) and end (red)\n' ...
    'Orbit markers: %d evenly spaced samples (requested)\n' ...
    'B-field model (simulation): %s\n' ...
    '\nStart snapshot:\n%s\n' ...
    '\nEnd snapshot:\n%s\n' ...
    '\nKey point: MTQ torque is always perpendicular to B.\n' ...
    'We visualize this as a plane normal to B at the satellite position.\n' ...
    ], numel(idx_mark), bfield_label(P), txt0, txtf);
annotation(fig,'textbox',[0.02 0.02 0.57 0.26], 'String', txt, ...
    'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], ...
    'Interpreter','none');

end

function [tau_dotB, ang_deg, txt] = plot_vectors_at_k(ax, P, S, k, L, showLegend)
r = S.r_eci(k,:).';

qk = S.q_ib(k,:);
C_bi = quat2dcm(qk);
C_ib = C_bi.'; % body->inertial

B_eci = S.B_eci(k,:).';
Bhat = B_eci / max(norm(B_eci), 1e-30);

m_body = S.m(k,:).';
m_eci = C_ib * m_body;

h_body = S.h_w(k,:).';
h_eci = C_ib * h_body;

tau_mtq_body = S.tau_mtq_body(k,:).';
tau_mtq_eci = C_ib * tau_mtq_body;

tau_des_body = -P.KH * h_body;
tau_des_eci = C_ib * tau_des_body;

tau_par_eci = (dot(tau_des_eci, Bhat) * Bhat);
tau_perp_eci = tau_des_eci - tau_par_eci;

tau_dotB = dot(tau_mtq_eci, Bhat);
ang_deg = acosd(max(-1, min(1, dot(tau_mtq_eci, Bhat) / max(norm(tau_mtq_eci), 1e-30))));

hv = 'off';
if showLegend; hv = 'on'; end

% Plane perpendicular to B at satellite
draw_perp_plane(ax, r, Bhat, 0.14*P.earth.Re_m, showLegend);

% Vectors (normalized directions)
plot_vec(ax, r, B_eci, L, [0.0 0.45 0.74], 'B', '-', hv);
plot_vec(ax, r, h_eci, L, [0.35 0.35 0.35], 'h_w (wheel momentum)', '-', hv);
plot_vec(ax, r, m_eci, L, [0.85 0.33 0.10], 'm (dipole)', '-', hv);
plot_vec(ax, r, tau_mtq_eci, L, [0.47 0.67 0.19], '\tau_{mtq}=m×B', '-', hv);

% Offset the decomposition vectors so tau_mtq (often == tau_des,⊥B) does not
% visually hide the perpendicular/parallel components.
satSizeRe = 0.02;
if isfield(P,'viz') && isfield(P.viz,'sat_size_Re')
    satSizeRe = double(P.viz.sat_size_Re);
end
satOffset = 1.6 * (satSizeRe * P.earth.Re_m);
ex = C_ib(:,1);
ey = C_ib(:,2);
r_perp = r + satOffset*ex;
r_par  = r + satOffset*ey;
plot_vec(ax, r_perp, tau_perp_eci, L, [0 0 0], '\tau_{des,\perp B}', '--', hv);
plot_vec(ax, r_par,  tau_par_eci,  L, [0.7 0 0.7], '\tau_{des,\parallel B}', ':', hv);

% Body axes triad
draw_body_axes(ax, r, C_ib, 0.16*P.earth.Re_m, hv);

tag = "start";
if ~showLegend; tag = "end"; end
text(ax, r(1), r(2), r(3), ['  ' tag], 'Color','k', 'FontWeight','bold');

txt = sprintf([ ...
    't = %.1f s\n' ...
    '||h_w|| = %.3e [N·m·s]\n' ...
    '||B|| = %.1f [µT]\n' ...
    '||m|| = %.3f [A·m^2]\n' ...
    '||tau_{mtq}|| = %.2f [µN·m]\n' ...
    'tau_{mtq}·bhat = %.2e (ideal 0)\n' ...
    'angle(tau_{mtq},B) = %.2f deg (ideal 90)\n' ...
    ], S.t(k), norm(h_body), 1e6*norm(B_eci), norm(m_body), 1e6*norm(tau_mtq_body), tau_dotB, ang_deg);
end

function plot_vec(ax, r, v, L, color, name, lineStyle, hv)
[u, ok] = unit_vec(v);
if ~ok; u = [0;0;0]; end
h = quiver3(ax, r(1), r(2), r(3), L*u(1), L*u(2), L*u(3), 0, ...
    'Color',color, 'LineWidth', 1.25, 'MaxHeadSize', 0.7, ...
    'LineStyle',lineStyle, 'HandleVisibility', hv, 'DisplayName', name);
try, h.Clipping = 'off'; end
end

function draw_perp_plane(ax, r, n_hat, halfSize, showLegend)
% Draw a translucent plane patch with normal n_hat centered at r.
% The plane indicates the "allowed torque plane" (⊥B).

if nargin < 5; showLegend = false; end
hv = 'off';
if showLegend; hv = 'on'; end

% Choose spanning vectors
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

patch(ax, X, Y, Z, [0.6 0.6 0.6], ...
    'FaceAlpha', 0.12, 'EdgeColor',[0.5 0.5 0.5], 'LineStyle','--', ...
    'HandleVisibility', hv, 'DisplayName','plane ⟂B (achievable τ)');
end

function draw_body_axes(ax, r, C_ib, L, hv)
% Draw body axes triad in ECI (x=red,y=green,z=blue)

ex = C_ib(:,1); ey = C_ib(:,2); ez = C_ib(:,3);
hx = quiver3(ax, r(1), r(2), r(3), L*ex(1), L*ex(2), L*ex(3), 0, ...
    'Color',[0.85 0.33 0.10], 'LineWidth', 1.0, 'MaxHeadSize', 0.5, ...
    'HandleVisibility', hv, 'DisplayName','body x');
hy = quiver3(ax, r(1), r(2), r(3), L*ey(1), L*ey(2), L*ey(3), 0, ...
    'Color',[0.25 0.6 0.2], 'LineWidth', 1.0, 'MaxHeadSize', 0.5, ...
    'HandleVisibility', hv, 'DisplayName','body y');
hz = quiver3(ax, r(1), r(2), r(3), L*ez(1), L*ez(2), L*ez(3), 0, ...
    'Color',[0.1 0.3 0.8], 'LineWidth', 1.0, 'MaxHeadSize', 0.5, ...
    'HandleVisibility', hv, 'DisplayName','body z');
try, hx.Clipping = 'off'; hy.Clipping = 'off'; hz.Clipping = 'off'; end
end

function draw_sat_cube(ax, r_eci, q_ib, size_m, color, alpha, name)
% Draw a simple cube representing the spacecraft, oriented by q_ib.

C_bi = quat2dcm(q_ib);
C_ib = C_bi.';

s = 0.5 * size_m;
Vb = [ ...
    -s -s -s;
     s -s -s;
     s  s -s;
    -s  s -s;
    -s -s  s;
     s -s  s;
     s  s  s;
    -s  s  s ].';

Vi = (C_ib * Vb) + r_eci;
V = Vi.';

F = [ ...
    1 2 3 4; % bottom
    5 6 7 8; % top
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8 ];

h = patch(ax, 'Vertices', V, 'Faces', F, ...
    'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor',[1 1 1]*0.9, ...
    'LineWidth', 0.4, 'DisplayName', name);
try, h.Clipping = 'off'; end
end

function set_camtarget_to_cursor(ax)
try
    cp = ax.CurrentPoint;
    ax.CameraTarget = cp(1,1:3);
catch
end
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
% Analytic dipole field lines (for visualization context).

Ls = [1.6 2.0 3.0 4.5];
phis = deg2rad([0 60 120 180 240 300]);
theta = linspace(0.15, pi-0.15, 400); % colatitude

tilt = 0;
if isfield(P,'earth') && isfield(P.earth,'dipole_tilt_deg')
    tilt = deg2rad(P.earth.dipole_tilt_deg);
end
Ry = [cos(tilt) 0 sin(tilt); 0 1 0; -sin(tilt) 0 cos(tilt)]; % tilt about +Y (X-Z plane)

for L = Ls
    rr = L * (sin(theta).^2);
    for phi = phis
        x = rr .* sin(theta) * cos(phi);
        y = rr .* sin(theta) * sin(phi);
        z = rr .* cos(theta);
        pts = Ry * [Re*x; Re*y; Re*z];
        plot3(ax, pts(1,:), pts(2,:), pts(3,:), '-', 'Color',[0.2 0.7 0.9], 'LineWidth', 0.6, 'HandleVisibility','off');

        % Direction arrows
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
