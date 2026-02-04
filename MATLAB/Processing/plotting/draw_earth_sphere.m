function h = draw_earth_sphere(ax, Re_m, P)
%DRAW_EARTH_SPHERE Draw a nicer-looking Earth sphere (optionally textured).
%
% Uses built-in demo data (topo/topomap1) if available, otherwise falls back
% to a solid-color sphere. Lighting is enabled for 3D depth.

if nargin < 3
    P = struct();
end

% Resolution (keep moderate for interaction performance).
n = 100;
if isfield(P,'viz') && isfield(P.viz,'earth_res')
    n = max(30, round(P.viz.earth_res));
end

[xe, ye, ze] = sphere(n);

% Default appearance
alpha = 0.95;
if isfield(P,'viz') && isfield(P.viz,'earth_alpha')
    alpha = max(0.05, min(1.0, double(P.viz.earth_alpha)));
end

% Try texture from built-in topo dataset
use_tex = true;
if isfield(P,'viz') && isfield(P.viz,'earth_texture')
    use_tex = logical(P.viz.earth_texture);
end

if use_tex
    try
        load topo topo topomap1 %#ok<LOAD>
        topo = double(topo);
        % Match topo to sphere grid size (no extra toolboxes).
        ny = size(xe,1);
        nx = size(xe,2);
        [Xq, Yq] = meshgrid(linspace(1, size(topo,2), nx), linspace(1, size(topo,1), ny));
        topo_rs = interp2(topo, Xq, Yq, 'linear');

        h = surf(ax, Re_m*xe, Re_m*ye, Re_m*ze, topo_rs, ...
            'FaceColor','texturemap', 'EdgeColor','none', 'FaceAlpha', alpha, ...
            'HandleVisibility','off');
        try
            colormap(ax, topomap1);
        catch
            colormap(ax, parula(256));
        end
    catch
        h = [];
    end
else
    h = [];
end

if isempty(h)
    h = surf(ax, Re_m*xe, Re_m*ye, Re_m*ze, ...
        'FaceColor',[0.15 0.35 0.7], 'EdgeColor','none', 'FaceAlpha', alpha, ...
        'HandleVisibility','off');
end

% Lighting defaults
try
    set(h, 'FaceLighting','gouraud', 'SpecularStrength',0.25, 'DiffuseStrength',0.85, 'AmbientStrength',0.25);
    shading(ax,'interp');
catch
end
end

