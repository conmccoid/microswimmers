%% plot_fd_rotlet_3d.m
%
% Visualises a single force-doublet PUSHER at the origin together with
% the contribution of a co-located rotlet.
%
% Layout
%   - 3-D figure with a pcolor slice (z = 0 plane) showing |u_fd| (force
%     doublet speed magnitude, colour-mapped)
%   - Streamlines in the z = 0 plane showing the ROTLET velocity field,
%     drawn over the colour slice
%   - A red arrow indicating the swimmer orientation (beta = [1 0 0])
%   - A semi-transparent sphere at the origin representing the swimmer body
%
% All kernels are the Cortez/Ainley regularised (blob) Stokeslet-derivative
% kernels — no external helper needed.
%
% Parameters you may want to tweak
%   d    : regularisation length (blob radius)
%   q0   : stresslet  / force-doublet strength (pusher => q0 > 0)
%   b0   : source-dipole strength
%   r0   : rotlet strength
%   Cr   : rotlet asymmetry parameter
%   L    : half-width of the plotting domain (in units of d)
%   n    : grid resolution (n x n points per slice)

clear; close all; clc

%% ── Parameters ──────────────────────────────────────────────────────────
d   = 0.2;          % blob radius
ppfac = +1.0;       % +1 pusher, -1 puller

q0  =  4*pi*ppfac;  % stresslet strength  (matches fd_script convention)
U   =  1.0;         % self-swimming speed
b0  = -2*pi*d^3*U;  % source-dipole strength

r0  =  6*pi;        % rotlet strength
Cr  = -1.0;         % rotlet asymmetry

beta = [1 0 0];     % swimmer orientation (unit vector)
xs   = [0 0 0];     % swimmer position

L  = 1.2;           % half-domain (physical units)
n  = 120;           % grid points per axis (increase for sharper images)

%% ── Build 2-D slice grid (z = 0) ────────────────────────────────────────
x_grid = linspace(-L, L, n);
[Xsl, Ysl] = meshgrid(x_grid, x_grid);
Zsl = zeros(size(Xsl));

% Flatten for vectorised kernel evaluation
pts = [Xsl(:) - xs(1), Ysl(:) - xs(2), Zsl(:) - xs(3)];  % (n^2, 3)
Np  = size(pts, 1);


%% ── Evaluate fields on the slice ────────────────────────────────────────
U_fd  = fd_kernel(pts,  beta, d, q0, b0);
U_rot = rotlet_kernel(pts, beta, d, r0, Cr);

% Speed magnitude of force-doublet field (for colour)
spd_fd  = reshape(vecnorm(U_fd,  2, 2), n, n);

% Rotlet in-plane components for streamlines
Ux_rot = reshape(U_rot(:,1), n, n);
Uy_rot = reshape(U_rot(:,2), n, n);

%% ── Also compute on 3-D volume for slice3 surface ───────────────────────
% (same z=0 plane — we embed it in a thin volume for slice() compatibility)
% We use the 2-D results directly via surf/pcolor in 3-D axes.

%% ── Figure ───────────────────────────────────────────────────────────────
fig = figure('Color','k','Position',[100 100 1100 860]);

ax = axes('Parent', fig);
set(ax, 'Color','k', 'GridColor','w', 'XColor','w', 'YColor','w', 'ZColor','w', ...
    'FontSize', 22, 'TickLabelInterpreter', 'latex', 'LineWidth', 1.2)
hold on; grid on; box on;

%── (1) pcolor slice ──────────────────────────────────────────────────────
hsurf = surf(Xsl, Ysl, Zsl, spd_fd, 'EdgeColor','none', 'FaceAlpha', 1.0);

% Cap colourscale so the near-field singularity does not swamp the palette
clim_max = prctile(spd_fd(:), 97);
clim([0, clim_max])

colormap(ax, parula)
cb = colorbar('Color','w','TickLabelInterpreter','latex','FontSize',18);
cb.Label.String = '$|\mathbf{u}_{FD}|$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 22;

%── (2) Rotlet streamlines in z = 0 plane ────────────────────────────────
% Seed points along a circle around the origin
n_seeds = 32;
theta_s = linspace(0, 2*pi, n_seeds+1); theta_s(end) = [];
r_seed  = 0.55;
sx = r_seed * cos(theta_s);
sy = r_seed * sin(theta_s);

% streamline() needs 3-D arrays; expand the 2-D slice into a 1-layer volume
Xs3 = reshape(Xsl,  [n n 1]);
Ys3 = reshape(Ysl,  [n n 1]);
Zs3 = zeros([n n 1]);
Ux3 = reshape(Ux_rot, [n n 1]);
Uy3 = reshape(Uy_rot, [n n 1]);
Uz3 = zeros([n n 1]);

hSL = streamline(Xs3, Ys3, Zs3, Ux3, Uy3, Uz3, sx(:), sy(:), zeros(size(sx(:))));
set(hSL, 'Color', [1 1 1 0.75], 'LineWidth', 1.6)

%── (3) Swimmer body sphere ───────────────────────────────────────────────
[sphX, sphY, sphZ] = sphere(40);
r_body = 0.85*d;
hSph = surf(r_body*sphX, r_body*sphY, r_body*sphZ, ...
            'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'none', 'FaceAlpha', 0.85);

%── (4) Orientation arrow ─────────────────────────────────────────────────
arr_len = 3.5*d;
quiver3(0, 0, 0, arr_len*beta(1), arr_len*beta(2), arr_len*beta(3), ...
        0, 'r', 'LineWidth', 3.5, 'MaxHeadSize', 1.2)

%── (5) Labels & limits ───────────────────────────────────────────────────
xlabel('$x$', 'Interpreter','latex','Color','w')
ylabel('$y$', 'Interpreter','latex','Color','w')
zlabel('$z$', 'Interpreter','latex','Color','w')

if ppfac > 0
    type_str = 'Pusher';
else
    type_str = 'Puller';
end
title(sprintf('%s: force-doublet $|\\mathbf{u}_{FD}|$ with rotlet streamlines', type_str), ...
      'Interpreter','latex','Color','w','FontSize',20)

xlim([-L L]); ylim([-L L]); zlim([-0.5 0.5])
view([-30 35])
axis equal tight
zlim([-0.3 0.3])
daspect([1 1 0.35])   % flatten slightly so the 2-D slice is prominent

%── (6) Optional: add a faint vertical plane (x = 0) to give 3-D depth ───
Xv = zeros(2,2);
Yv = [-L  L; -L  L];
Zv = [-0.3 -0.3; 0.3  0.3];
surf(Xv, Yv, Zv, 'FaceAlpha', 0.08, 'FaceColor', 'w', 'EdgeColor', 'none')

hold off

%% ── Export ───────────────────────────────────────────────────────────────
set(gcf, 'PaperUnits','points', 'PaperPositionMode','manual', ...
         'PaperSize',[1200,1000], 'PaperPosition',[0 0 1100 860])

print('fd_pusher_rotlet_3d', '-dpng', '-r300')
print('fd_pusher_rotlet_3d', '-depsc', '-r300')
fprintf('Saved fd_pusher_rotlet_3d.png / .eps\n')

%% ── Regularised kernel helpers (inline, no external files needed) ───────

    function u = fd_kernel(pts, beta, d, q0, b0)
        % Vectorised force-doublet velocity at pts due to one swimmer at origin
        r   = vecnorm(pts, 2, 2);
        Np  = size(pts, 1);
        bdx = pts * beta';                      % beta · x  (N,1)
        R   = sqrt(r.^2 + d^2);
        R3  = R.^3;   R5 = R.^5;
        F1  = (1/8/pi) * (-pts./R3 + 3*(bdx.^2).*pts./R5 ...
                          + 3*bdx.*d^2.*(ones(Np,1)*beta)./R5);
        F2  = (1/8/pi) * (-3*bdx.*pts./R5 + (R.^2 - 3*d^2).*(ones(Np,1)*beta)./R5);
        u   = q0.*F1 + b0.*F2;
        % zero out the self-point (r ~ 0)
        u(r < 1e-12, :) = 0;
    end

    function u = rotlet_kernel(pts, beta, d, r0, Cr)
        % Vectorised rotlet velocity at pts due to one swimmer at origin
        r    = vecnorm(pts, 2, 2);
        Np  = size(pts, 1);
        % cross(beta, x) for each row
        rxx  = cross(ones(Np,1)*beta, pts, 2);
        rdx  = pts * beta';                     % beta · x  (N,1)
        R    = sqrt(r.^2 + d^2);
        R5   = R.^5;   R7 = R.^7;
        F1   = (1/8/pi) .* rdx .* rxx .* (6./R5 + 15*d^2./R7);
        F2   = (1/8/pi) .* Cr  .* rxx .* 15*d^3 ./ R7;
        u    = r0.*(F1 - F2);
        u(r < 1e-12, :) = 0;
    end
