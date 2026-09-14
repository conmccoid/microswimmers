
%% This code is for plotting the resulting flow field for one force doublet using a 7/2 blob


n=101;
L=1.0;
nskip=4.0;
x_grid=linspace(-L,L,n);
[X,Y,Z]=meshgrid(x_grid,x_grid,x_grid);
vX=0*Z;
vY=0*Y;
vZ=0*X;
vx=zeros(n,1);
vy=zeros(n,1);
vz=zeros(n,1);  

ppfac=+1.0;
d = 0.2;
 
q0 =  1.0*ppfac;
p0 = 1*q0;
U0=p0*(4*pi*d);
p0=(4*pi*d);
d3=d*d*d;
d5=d*d*d3;
d7=d*d*d5;


xs=zeros(3,1);
beta=[1 0 0];
n3=n^3;
for j=1:n
    for k=1:n
        for el=1:n
            
            
            dx=X(j,k,el)-xs(1);
            dy=Y(j,k,el)-xs(2);
            dz=Z(j,k,el)-xs(3);
            
            betadotx=beta(1)*dx+beta(2)*dy+beta(3)*dz;
            betadotx2=betadotx*betadotx;
            z1=sqrt(dx*dx+dy*dy+dz*dz)/d;
            z2 = z1*z1;
            z3 = z1*z2;
            z4 = z1*z3;
            z5 = z1*z4;
            z6 = z1*z5;
            z7 = z1*z6;
            % zexp=exp(-z1);
            r=sqrt(dx*dx+dy*dy+dz*dz);
            [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d);
            
            vX(j,k,el)=-q0*beta(1)*(J1)*betadotx- q0*dx*(H2+betadotx2*J0) -p0*(beta(1)*J2+betadotx*dx*J3);
            vY(j,k,el)=-q0*beta(2)*(J1)*betadotx- q0*dy*(H2+betadotx2*J0) -p0*(beta(2)*J2+betadotx*dy*J3);
            vZ(j,k,el)=-q0*beta(3)*(J1)*betadotx- q0*dz*(H2+betadotx2*J0) -p0*(beta(3)*J2+betadotx*dz*J3);
        end
        
    end
end
vMag=sqrt(vX.^2+vY.^2+vZ.^2);

figure(1)
% quiver(X(:,:,ceil(n/2)),Y(:,:,ceil(n/2)),vX(:,:,ceil(n/2)),vY(:,:,ceil(n/2)),'LineWidth',1.0,'MaxHeadSize',5.0,'AutoScale', 'on')
% quiver(X(1:nskip:n,1:nskip:n,ceil(n/2)),Y(1:nskip:n,1:nskip:n,ceil(n/2)),vX(1:nskip:n,1:nskip:n,ceil(n/2)),vY(1:nskip:n,1:nskip:n,ceil(n/2)),'LineWidth',1.0,'MaxHeadSize',1.0,'AutoScale', 'on','Color','k')

% contour(X(:,:,ceil(n/2)),Y(:,:,ceil(n/2)),vMag(:,:,ceil(n/2)),30,'LineWidth', 1.50)
pcolor(X(:,:,ceil(n/2)),Y(:,:,ceil(n/2)),vMag(:,:,ceil(n/2)))
shading interp;
% [sX,sY,sZ]=meshgrid([-L L],x_grid(1:10:end),0);
hold on
hstream=streamslice(X,Y,Z,vX,vY,vZ,[],[],0,2,'arrows','nearest');
set(hstream,'color','w','LineWidth',2)
% quiver(0,0,.35,0,'LineWidth',3.0,'MaxHeadSize',5.0,'Color','r')
Arrow0 = annotation('arrow') ;
Arrow0.Parent = gca; 
Arrow0.Position = [-0.05, 0, .1, 0] ; %puller
Arrow0.Color='red';
Arrow0.LineWidth=3.0;
Arrow0.HeadWidth=30.0;
hold off
contourcbar

set(gca, 'FontSize', 30,'TickLabelInterpreter','latex')
% title('Pusher, \epsilon=0.1, C=-0.1')
% contourcbar
caxis([0 1])
set(gcf,'paperunits','points')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf,'PaperSize',[1200,1000])
set(gcf,'PaperPosition',[0,0,1000,800])
% axis([.1 100 0. 5.])
% print(['pusher_contours_pcolor'],'-dpng','-r300');
% print(['pusher_contours_pcolor'],'-depsc','-r300');
print(['72blob_contours_pcolor'],'-dpng','-r300');
print(['72blob_contours_pcolor'],'-depsc','-r300');
