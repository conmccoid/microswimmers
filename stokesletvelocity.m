function u_s = stokesletvelocity(x_fd,x_s,f_s,d)
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

u_s=zeros(size(x_fd));
N=size(x_fd,1);
M=size(x_s,1);
for i=1:N
    r=vecnorm(x_fd(i,:)-x_s,2,2);

    x=x_fd(i,:)-x_s;
    
    % dx=x_fd(i,1)-x_s(:,1);
    % dy=x_fd(i,2)-x_s(:,2);
    % dz=x_fd(i,3)-x_s(:,3);
    % fdotx=f_s(:,1).*dx+f_s(:,2).*dy+f_s(:,3).*dz;
    % 
    fdotx=dot(f_s,x,2);

    R=(d.^2+r.^2).^(1/2);
    H1=(1/8).*pi.^(-1).*(R.^(-1)+d.^(2).*R.^(-3));
    H2=(1/8).*pi.^(-1).*R.^(-3);
    Fact1=f_s.*H1;
    Fact2=fdotx.*H2.*x;
    u_s(i,:)=sum(Fact1+Fact2,1)';
    % [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d);
    % for j=1:M
    %     u_s(i,1)=u_s(i,1)+f_s(j,1)*H1(j)+ fdotx(j)*dx(j)*H2(j);
    %     u_s(i,2)=u_s(i,2)+f_s(j,2)*H1(j)+ fdotx(j)*dy(j)*H2(j);
    %     u_s(i,3)=u_s(i,3)+f_s(j,3)*H1(j)+ fdotx(j)*dz(j)*H2(j);
    % 
    % end
end