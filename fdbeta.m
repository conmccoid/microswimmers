function u_betafd = fdbeta(x_fd,beta_fd,d,q0,b0)
% The goal of this function is to compute the force doublet velocities
% for each force doublet position. 
% This is NOT for grid evaluation
% Since this is Stokes regime, each force doublet has a) an intrinsic
% swimming speed (based on del/C/q0,p0 and current orientation beta);
% and b) contributions from the other force doublets
% Note: We do this by first creating a copy of each location
% INPUTS
% x_s = location of point forces
% beta = orientation of force doublets
% d = regularization parameters
% f_s = forces
% OUTPUT
% u_s = velocity of force doublet due to stokeslets
% 

x_sd=x_fd;
beta_sd=beta_fd;

u_betafd=zeros(size(x_fd));
N=size(x_fd,1);
% M=size(x_sd,1);
for i=1:N
    r=vecnorm(x_fd(i,:)-x_sd,2,2);
    
    x=x_fd(i,:)-x_sd;
    
    bxx=cross(ones(N,1)*beta_fd(i,:),x,2);
    bxxxb=cross(bxx,ones(N,1)*beta_fd(i,:),2);
    bxr=cross(ones(N,1)*beta_fd(i,:),beta_sd,2);
    bxrxb=cross(bxr,ones(N,1)*beta_fd(i,:),2);
    
    bdotr=dot(ones(N,1)*beta_fd(i,:),beta_sd,2);
    rdotx=dot(beta_sd,x,2);
    rdotx2=rdotx.^2;
    bdotx=dot(ones(N,1)*beta_fd(i,:),x,2);
    
    
    R=(d.^2+r.^2).^(1/2);
    R2=R.*R;
    R3=R.*R2;
    R5=R3.*R2;
    R7=R5.*R2;
    
    J1=6./R5.*bdotr.*rdotx+3.*(1./R5-5.*rdotx2./R7).*bdotx;
    J2=3*d.^(2).*(bdotr./R5-5.*rdotx.*bdotx./R7);
    J3=-3./R5.*(bdotr)+15./R7.*rdotx.*bdotx;
    J4=(-3./R5+d.^(2).*15./R7).*bdotx;
    Fact1=(1/8).*pi.^(-1).*bxxxb.*(-q0.*J1+b0.*J3);
    Fact2=(1/8).*pi.^(-1).*bxrxb.*(-q0.*J2+b0.*J4);
    u_betafd(i,:)=sum(Fact1+Fact2,1)';
end
