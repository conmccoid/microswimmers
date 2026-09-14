function u_fd = fdvelocity(x_fd,beta_fd,d,q0,b0)
% The goal of this function is to compute the force doublet velocities
% for each force doublet position. 
% This is NOT for grid evaluation
% Since this is Stokes regime, each force doublet has a) an intrinsic
% swimming speed (based on del/C/q0,p0 and current orientation beta);
% and b) contributions from the other force doublets
% Note: We do this by first creating a copy of each location
% INPUTS
% x_fs = location of force doublets
% beta = orientation of force doublets
% del = regularization parameters
% asymmetry parameters 
% q0, b0 = dipole and stresslet strength
% OUTPUT
% u_fd = velocity of force doublet due to force doublets


u_fd=zeros(size(x_fd));
x_cp=x_fd;
N=size(x_fd,1);

for i=1:N
    r=vecnorm(x_fd(i,:)-x_cp,2,2);
    x=x_fd(i,:)-x_cp;
    % dx=x_fd(i,1)-x_cp(:,1);
    % dy=x_fd(i,2)-x_cp(:,2);
    % dz=x_fd(i,3)-x_cp(:,3);

    % betadotx=beta(:,1).*dx+beta(:,2).*dy+beta(:,3).*dz;
    betadotx=dot(ones(N,1)*beta_fd(i,:),x,2);
    betadotx2=betadotx.*betadotx;
    R=sqrt(r.^2+d.^2);
    R2=R.*R;
    R3=R2.*R;
    R5=R3.*R2;
    Fact1=(1/8).*pi^(-1).*(-x./R3+3.*betadotx2.*x./R5+3.*betadotx.*d.^2.*ones(N,1)*beta_fd(i,:)./R5);
    Fact2=(1/8).*pi^(-1).*(-3.*betadotx.*x./R5+(R2-3.*d.^2).*ones(N,1)*beta_fd(i,:)./R5);
    
    u_fd(i,:)=sum(q0.*Fact1+b0.*Fact2,1)';    
    % Fact1=(1/8).*pi.^(-1)*-1.*R^(-3)+3.*betadotx2.*R^(-5);
    % [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d);
    % for j=1:N
    %     u_fd(i,1)=u_fd(i,1)+q0*beta(j,1)*(J1(j))*betadotx(j)+ q0*dx(j)*(H2(j)+betadotx2(j)*J0(j)) +b0*(beta(i,1)*J2(j)+betadotx(j)*dx(j)*J3(j));
    %     u_fd(i,2)=u_fd(i,2)+q0*beta(j,2)*(J1(j))*betadotx(j)+ q0*dy(j)*(H2(j)+betadotx2(j)*J0(j)) +b0*(beta(i,2)*J2(j)+betadotx(j)*dy(j)*J3(j));
    %     u_fd(i,3)=u_fd(i,3)+q0*beta(j,3)*(J1(j))*betadotx(j)+ q0*dz(j)*(H2(j)+betadotx2(j)*J0(j)) +b0*(beta(i,3)*J2(j)+betadotx(j)*dz(j)*J3(j));
    % end
end