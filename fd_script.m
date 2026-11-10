%% force doublet script
clear all
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
    end
end

profile viewer
profile off