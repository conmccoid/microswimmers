function f_s= stericforce(x_fd, kappa, tol)
% INPUTS
% x_fd = locations of force doublets (nv,3)
% kappa = steric force strength i.e. kappa(x_fd1-x_fd2)
% tol = tolerance distance when a force doublet is close to another 
% OUTPUTS
% f_s = steric forces (nv,3)

f_s=zeros(size(x_fd));

N=size(x_fd,1);

for i=1:N
    x=x_fd(i,:)-x_fd;
    x_dist=vecnorm(x_fd(i,:)-x_fd,2,2);
    [tol_idx] = find(x_dist<tol);
    tol_idx(tol_idx==i)=[];
    f_s(i,:)=sum(kappa.*(tol-x_dist(tol_idx)).*x(tol_idx,:)./x_dist(tol_idx),1);
    % f_s(i,:)=f_s(tol_idx,:)+kappa.*(tol-x_dist(tol_idx)).*x(tol_idx,:)./x_dist(tol_idx);
end