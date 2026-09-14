function write_pvd(j_max,ppfac)

if ppfac==1
    f_type='pusher';
elseif ppfac==-1
    f_type='puller';
else
    f_type='particle';
end

fname = sprintf([f_type,'_%05d.vtk'], j);
fid   = fopen(fname, 'w');
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="Collection" version="0.1">\n');
fprintf(fid, '  <Collection>\n');
for k = 1:j_max
    fprintf(fid, '    <DataSet timestep="%d" file="particles_%05d.vtk"/>\n', k, k);
end
fprintf(fid, '  </Collection>\n');
fprintf(fid, '</VTKFile>\n');
fclose(fid);
end