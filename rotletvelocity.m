function u_r = rotletvelocity(x_fd,beta_fd,d,r0,Cr)



x_sd=x_fd;
beta_sd=beta_fd;

u_r=zeros(size(x_fd));
N=size(x_fd,1);

for i=1:N
    r=vecnorm(x_fd(i,:)-x_sd,2,2);

    x=x_fd(i,:)-x_sd;

    rxx=cross(beta_sd,x,2);    
    rdotx=dot(beta_sd,x,2);

    R=(d.^2+r.^2).^(1/2);
    R5=R.^(5);
    R7=R.^(7);

    Fact1=(1/8).*pi.^(-1).*rdotx.*rxx.*(6./R5+15.*d.^2./R7);
    Fact2=(1/8).*pi.^(-1).*Cr.*rxx.*15.*d.^3./R7;
    u_r(i,:)=sum(r0.*Fact1-r0.*Fact2,1)';
end


end