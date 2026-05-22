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
    
    % dx=x_fd(i,1)-x_s(:,1);
    % dy=x_fd(i,2)-x_s(:,2);
    % dz=x_fd(i,3)-x_s(:,3);
    % fdotx=f_s(:,1).*dx+f_s(:,2).*dy+f_s(:,3).*dz;
    % 
    bxx=cross(beta_fd,x,2);
    bxxxb=cross(bxx,beta_fd,2);
    bxr=cross(beta_fd,beta_sd,2);
    bxrxb=cross(bxr,beta_fd,2);
    
    bdotr=dot(beta_fd,beta_sd,2);
    rdotx=dot(beta_sd,x,2);
    rdotx2=rdotx.^2;
    bdotx=dot(beta_fd,x,2);
    
    
    R=(d.^2+r.^2).^(1/2);
    J1=6.*R.^(-5).*bdotr.*rdotx+3.*(R.^(-5)-5.*rdotx2.*R.^(-7)).*bdotx;
    J2=3*d.^(2).*(bdotr.*R.^(-5)-5.*rdotx.*bdotx.*R.^(-7));
    J3=-3.*R.^(-5).*(bdotr)+15.*R.^(-7).*rdotx.*bdotx;
    J4=(-3.*R.^(-5)+d.^(2).*15.*R.^(-7)).*bdotx;
    Fact1=(1/8).*pi.^(-1).*bxxxb.*(q0.*J1+b0.*J3);
    Fact2=(1/8).*pi.^(-1).*bxrxb.*(q0.*J2+b0.*J4);
    u_betafd(i,:)=sum(Fact1+Fact2,1)';
    % for j=1:M
    %     betadotx=beta(i,1)*dx(j)+
    %     u_betas(i,1)=u_betas(i,1)+f_s(j,1)*H1(j)+ fdotx(j)*dx(j)*H2(j);
    %     u_betas(i,2)=u_betas(i,2)+f_s(j,2)*H1(j)+ fdotx(j)*dy(j)*H2(j);
    %     u_betas(i,3)=u_betas(i,3)+f_s(j,3)*H1(j)+ fdotx(j)*dz(j)*H2(j);
    % 
    % end
end
