%% force doublet script
clear
profile on

d=1.0; %regularization parameter


% T=20.00;
T = 0.1;
dt=1e-3;
Nt=ceil(T/dt);
ppfac=-1.0; % pusher = 1, puller = -1.0


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
x_list=[ X(:) Y(:) Z(:)];


% All pointing the same way
% beta_fd=[ones(size(X(:))) zeros(size(X(:)))  zeros(size(X(:)))];
%
%Setting random orientation
theta=2*pi*rand(size(X(:)));
cosphi=rand(size(X(:)))*2-1;
phi=acos(cosphi);
beta_list=[ sin(phi).*cos(theta) sin(phi).*sin(theta) cos(phi) ];

x_store=zeros([size(x_list) ceil(Nt/nstore)]);
beta_store=zeros([size(x_list) ceil(Nt/nstore)]);
x_store(:,:,1)=x_list;
beta_store(:,:,1)=beta_list;
j=1;

x_store(:,:,j)    = x_list;
beta_store(:,:,j) = beta_list;
npts = size(x_list, 1);
f_s=zeros(size(x_list));
u_s=zeros(size(x_list));
u_fd=zeros(size(x_list));
u_r=zeros(size(x_list));
beta_s=zeros(size(x_list));
beta_fd=zeros(size(x_list));
beta_r=zeros(size(x_list));

write_vtk(j, x_list, beta_list, f_s, u_s, u_fd, u_r, beta_s, beta_fd, beta_r, ppfac);

for n=1:Nt
    %compute forces due to steric forces
    f_s=stericforce(x_list, kappa, tol);

    % update velocities and directions
    for i=1:npts
        [u_s(i,:), beta_s(i,:), u_fd(i,:), beta_fd(i,:), u_r(i,:), beta_r(i,:)] = updateSwimmer(i,x_list,beta_list, f_s, d, q0, b0, r0, Cr);
    end

    x_list=x_list+dt.*(u_s+u_fd+u_r);
    beta_list=beta_list+dt.*(beta_s+beta_fd+beta_r);
    beta_list=beta_list./(vecnorm(beta_list,2,2));
    
    if mod(n, nstore) == 0
        j = j + 1;
        x_store(:,:,j)    = x_list;
        beta_store(:,:,j) = beta_list;
        write_vtk(j, x_list, beta_list, f_s, u_s, u_fd, u_r, beta_s, beta_fd, beta_r, ppfac);
    end
end

profile viewer
profile off