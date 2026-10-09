function u_betas = stokesletbeta(x_fd,beta_fd,x_s,f_s,d)
% The goal of this function is to compute the force doublet velocities
% for each force doublet position. 
% This is NOT for grid evaluation
% Since this is Stokes regime, each force doublet has a) an intrinsic
% swimming speed (based on del/C/q0,p0 and current orientation beta);
% and b) contributions from the other force doublets
% Note: We do this by first creating a copy of each location
% INPUTS
% x_s = location of point forces (dim Mx3)
% beta = orientation of force doublets
% d = regularization parameters
% f_s = forces
% OUTPUT
% u_s = velocity of force doublet due to stokeslets
% 

u_betas=zeros(size(x_fd));
N=size(x_fd,1);
M=size(x_s,1); % number of point forces
for i=1:N
    r=vecnorm(x_fd(i,:)-x_s,2,2); % distance to point forces
    R=sqrt(r.^2+d.^2);
    R2=R.*R;
    R3=R2.*R;
    R5=R2.*R3;
    x=x_fd(i,:)-x_s;
    
    % dx=x_fd(i,1)-x_s(:,1);
    % dy=x_fd(i,2)-x_s(:,2);
    % dz=x_fd(i,3)-x_s(:,3);
    % fdotx=f_s(:,1).*dx+f_s(:,2).*dy+f_s(:,3).*dz;

    beta_i = beta_fd(i,:);
    bdotx = x*beta_i';
    bdotf = f_s*beta_i';
    bxxxb = norm(beta_i)*x - bdotx*beta_i;
    bxfxb = norm(beta_i)*f_s - bdotf*beta_i;
    % bxx=cross(ones(N,1)*beta_fd(i,:),x,2);
    % bxxxb=cross(bxx,ones(N,1)*beta_fd(i,:),2);
    % bxf=cross(ones(N,1)*beta_fd(i,:),f_s,2);
    % bxfxb=cross(bxf,ones(N,1)*beta_fd(i,:),2);
    
    % bdotf=dot(ones(N,1)*beta_fd(i,:),f_s,2);
    % bdotx=dot(ones(N,1)*beta_fd(i,:),x,2);
    fdotx=dot(f_s,x,2);
    
    Fact1=(1/8).*pi^(-1).*(bdotf./R3-3.*bdotx.*fdotx./R5).*bxxxb;
    Fact2=(1/8).*pi^(-1).*(1./R3+3.*d.^(2)./R5).*bdotx.*bxfxb;
    
    % [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d);

    u_betas(i,:)=sum(Fact1-Fact2,1)';
    % for j=1:M
    %     betadotx=beta(i,1)*dx(j)+
    %     u_betas(i,1)=u_betas(i,1)+f_s(j,1)*H1(j)+ fdotx(j)*dx(j)*H2(j);
    %     u_betas(i,2)=u_betas(i,2)+f_s(j,2)*H1(j)+ fdotx(j)*dy(j)*H2(j);
    %     u_betas(i,3)=u_betas(i,3)+f_s(j,3)*H1(j)+ fdotx(j)*dz(j)*H2(j);
    % 
    % end
end

%J3=(-3/8).*d.^2.*pi.^(-1).*R.^(-5);