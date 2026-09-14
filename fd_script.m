%% force doublet script
clear all

d=1.0; %regularization parameter


T=20.00;
dt=1e-3;
Nt=ceil(T/dt);
ppfac=0.0; % pusher = 1, puller = -1.0


q0=4*pi*ppfac; % stresslet strength
U=1.0; %self swimming velocity;
b0=-2*pi*d^3*U; %dipole
Cr=-1;
r0=6*pi; % rotlet
kappa=1.e3; %steric force strength
tol=1*d; %tolerance for when steric force is active

nstore=100;
t=0.0;

Nfd=8; %number of force doublet swimmers in a line
x_grid=linspace(-(Nfd-1)/2*2*tol,(Nfd-1)/2*2*tol,Nfd);

[X,Y,Z]=meshgrid(x_grid,x_grid, x_grid);
% [X,Y,Z]=meshgrid(x_grid,x_grid,x_grid);
x_fd=[ X(:) Y(:) Z(:)];


% All pointing the same way
% beta_fd=[ones(size(X(:))) zeros(size(X(:)))  zeros(size(X(:)))];
%
%Setting random orientation
theta=2*pi*rand(size(X(:)));
cosphi=rand(size(X(:)))*2-1;
phi=acos(cosphi);
beta_fd=[ sin(phi).*cos(theta) sin(phi).*sin(theta) cos(phi) ];



x_store=zeros([size(x_fd) ceil(Nt/nstore)]);
beta_store=zeros([size(x_fd) ceil(Nt/nstore)]);
x_store(:,:,1)=x_fd;
beta_store(:,:,1)=beta_fd;
j=1;

x_store(:,:,j)    = x_fd;
beta_store(:,:,j) = beta_fd;
npts = size(x_fd, 1);
f_s=zeros(size(x_fd));
u_s=zeros(size(x_fd));
u_fd=zeros(size(x_fd));
u_r=zeros(size(x_fd));
u_betas=zeros(size(x_fd));
u_betafd=zeros(size(x_fd));
u_betar=zeros(size(x_fd));

write_vtk(j, x_fd, beta_fd, f_s, u_s, u_fd, u_r, u_betas, u_betafd, u_betar, ppfac);

% if ppfac==1
%     f_type='pusher';
% elseif ppfac==-1
%     f_type='puller';
% else
%     f_type='particle';
% end
% 
% fname = sprintf([f_type,'_%05d.vtk'], j);
% fid   = fopen(fname, 'w');
% 
% 
% % --- Header ---
% fprintf(fid, '# vtk DataFile Version 3.0\n');
% fprintf(fid, 'Particles frame %d\n', j);
% fprintf(fid, 'ASCII\n');
% fprintf(fid, 'DATASET POLYDATA\n');
% 
% % --- Points ---
% fprintf(fid, 'POINTS %d float\n', npts);
% fprintf(fid, '%.6f %.6f %.6f\n', x_fd');   % transpose so fprintf walks rows
% 
% % --- Vertices (so ParaView renders them as glyphs) ---
% fprintf(fid, 'VERTICES %d %d\n', npts, 2*npts);
% for k = 0:npts-1
%     fprintf(fid, '1 %d\n', k);
% end
% 
% % --- Point data ---
% fprintf(fid, 'POINT_DATA %d\n', npts);
% 
% % orientation
% fprintf(fid, 'VECTORS orientation float\n');
% fprintf(fid, '%.6f %.6f %.6f\n', beta_fd');
% 
% % force
% fprintf(fid, 'VECTORS force float\n');
% fprintf(fid, '%.6f %.6f %.6f\n', f_s');
% 
% % velocity_stokeslet
% fprintf(fid, 'VECTORS velocity_stokeslet float\n');
% fprintf(fid, '%.6f %.6f %.6f\n', u_s');
% 
% % velocity_fd
% fprintf(fid, 'VECTORS velocity_fd float\n');
% fprintf(fid, '%.6f %.6f %.6f\n', u_fd');
% 
% fclose(fid);



for n=1:Nt
    %compute forces due to steric forces
    f_s=stericforce(x_fd, kappa, tol);
    % compute velocity and dbeta based off point forces 
    u_s = stokesletvelocity(x_fd,x_fd,f_s,d);
    u_betas = stokesletbeta(x_fd,beta_fd,x_fd,f_s,d);
    % compute velocity and dbeta based off force doublets
    u_fd = fdvelocity(x_fd,beta_fd,d,q0,b0);
    u_betafd = fdbeta(x_fd,beta_fd,d,q0,b0);
    % compute velocity and dbeta based off point forces 
    u_r = rotletvelocity(x_fd,beta_fd,d,r0,Cr);
    u_betar = rotletbeta(x_fd,beta_fd,d,r0,Cr);
    
    x_fd=x_fd+dt.*(u_s+u_fd+u_r);
    beta_fd=beta_fd+dt.*(u_betas+u_betafd+u_betar);
    beta_fd=beta_fd./(vecnorm(beta_fd,2,2));
    
    
    if mod(n, nstore) == 0
        j = j + 1;
        x_store(:,:,j)    = x_fd;
        beta_store(:,:,j) = beta_fd;
        write_vtk(j, x_fd, beta_fd, f_s, u_s, u_fd, u_r, u_betas, u_betafd, u_betar, ppfac);

    % 
    % fname = sprintf([f_type,'_%05d.vtk'], j);
    % fid   = fopen(fname, 'w');
    % 
    % npts = size(x_fd, 1);
    % 
    % % --- Header ---
    % fprintf(fid, '# vtk DataFile Version 3.0\n');
    % fprintf(fid, 'Particles frame %d\n', j);
    % fprintf(fid, 'ASCII\n');
    % fprintf(fid, 'DATASET POLYDATA\n');
    % 
    % % --- Points ---
    % fprintf(fid, 'POINTS %d float\n', npts);
    % fprintf(fid, '%.6f %.6f %.6f\n', x_fd');   % transpose so fprintf walks rows
    % 
    % % --- Vertices (so ParaView renders them as glyphs) ---
    % fprintf(fid, 'VERTICES %d %d\n', npts, 2*npts);
    % for k = 0:npts-1
    %     fprintf(fid, '1 %d\n', k);
    % end
    % 
    % % --- Point data ---
    % fprintf(fid, 'POINT_DATA %d\n', npts);
    % 
    % % orientation
    % fprintf(fid, 'VECTORS orientation float\n');
    % fprintf(fid, '%.6f %.6f %.6f\n', beta_fd');
    % 
    % % force
    % fprintf(fid, 'VECTORS force float\n');
    % fprintf(fid, '%.6f %.6f %.6f\n', f_s');
    % 
    % % velocity_stokeslet
    % fprintf(fid, 'VECTORS velocity_stokeslet float\n');
    % fprintf(fid, '%.6f %.6f %.6f\n', u_s');
    % 
    % % velocity_fd
    % fprintf(fid, 'VECTORS velocity_fd float\n');
    % fprintf(fid, '%.6f %.6f %.6f\n', u_fd');
    % 
    % fclose(fid);
    end
end


%%




%% Plotting

for k=1:j
    X = x_store(:,1,k);
    Y = x_store(:,2,k);
    Z = x_store(:,3,k);

    U = beta_store(:,1,k);
    V = beta_store(:,2,k);
    W = beta_store(:,3,k);

    clf
    figure(1)
    quiver3(X,Y,Z,U,V,W,'LineWidth',1.5)

    axis equal
    grid on

    xlabel('x')
    ylabel('y')
    zlabel('z')

    xlim([-1,7])
    ylim([-3,3])
    zlim([-3,3])
    title(['Frame ',num2str(k)])
    % pause
    drawnow


end