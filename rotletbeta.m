function u_betar = rotletbeta(x_fd,beta_fd,d,r0,Cr)


% The goal of this function is to compute the force doublet rotlet change
% in orientation
% for each force doublet position. 
% This is NOT for grid evaluation
% Note: We do this by first creating a copy of each location
% INPUTS
% x_fd = location of force doublets
% beta_fd = orientation of force doublets
% d = regularization parameter
% r0 = rotlet strength
% Cr = asymmetry parameter
% OUTPUT
% u_betar = velocity of force doublet due to stokeslets
% 

x_sd=x_fd;
beta_sd=beta_fd;

u_betar=zeros(size(x_fd));
N=size(x_fd,1);

for i=1:N
    r=vecnorm(x_fd(i,:)-x_sd,2,2);

    x=x_fd(i,:)-x_sd;

    rxx=cross(beta_sd,x,2);
    bxrxx=cross(ones(N,1)*beta_fd(i,:),rxx,2);
    bxrxxxb=cross(bxrxx,ones(N,1)*beta_fd(i,:),2);
    rxb=cross(beta_sd,ones(N,1)*beta_fd(i,:),2);

    bdotr=dot(ones(N,1)*beta_fd(i,:),beta_sd,2);
    bdotx=dot(ones(N,1)*beta_fd(i,:),x,2);
    rdotx=dot(beta_sd,x,2);

    R=(d.^2+r.^2).^(1/2);
    R2=R.*R;
    R3=R.*R2;
    R5=R3.*R2;
    R7=R5.*R2;
    R9=R7.*R2;

    Fact1=(1/8).*pi.^(-1).*bxrxxxb.*( (6./R5+15.*d.^2./R7).*bdotr -bdotx.*(rdotx.*(30./R7+105.*d.^2./R9)- Cr.*105.*d.^3./R9));
    Fact2=(1/8).*pi.^(-1).*rxb.*(rdotx.*(6./R5+15.*d.^2./R7) - Cr.*15.*d.^3./R7);
    u_betar(i,:)=sum(r0.*Fact1+r0.*Fact2,1)';
end
