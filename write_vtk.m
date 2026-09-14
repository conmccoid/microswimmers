function write_vtk(j, x_fd, beta_fd, f_s, u_s, u_fd, u_r, u_betas, u_betafd, u_betar, ppfac)

if ppfac==1
    f_type='pusher';
elseif ppfac==-1
    f_type='puller';
else
    f_type='particle';
end

vis_dir = fullfile('visualization', f_type);

if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end

fname = fullfile(vis_dir, sprintf('%s_%05d.vtk', f_type, j));
fid   = fopen(fname, 'w');
npts  = size(x_fd, 1);

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'Particles frame %d\n', j);
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');

fprintf(fid, 'POINTS %d float\n', npts);
fprintf(fid, '%.6f %.6f %.6f\n', x_fd');

fprintf(fid, 'VERTICES %d %d\n', npts, 2*npts);
for k = 0:npts-1
    fprintf(fid, '1 %d\n', k);
end

fprintf(fid, 'POINT_DATA %d\n', npts);

fields = {'orientation', beta_fd; ...
    'force',       f_s;     ...
    'velocity_stokeslet', u_s; ...
    'velocity_fd',        u_fd; ...
    'velocity_rotlet', u_r; ...
    'velocity_tot',  u_s+u_fd+u_r;...
    'dbeta_stokeslet', u_betas; ...
    'dbeta_fd',        u_betafd; ...
    'dbeta_rotlet', u_betar; ...
    'dbeta_tot',  u_betas+u_betafd+u_betar;...
    };

for f = 1:size(fields, 1)
    fprintf(fid, 'VECTORS %s float\n', fields{f,1});
    fprintf(fid, '%.6f %.6f %.6f\n',  fields{f,2}');
end

fclose(fid);
end