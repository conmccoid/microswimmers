%% force doublet script
clear all

d=0.2;
Nfd=11;
x_fd=[ zeros(Nfd,1) linspace(-0.5,0.5,Nfd)' zeros(Nfd,1)];

beta_fd=[ones(Nfd,1) zeros(Nfd,1)  zeros(Nfd,1)];

T=1.0;
dt=1e-4;
Nt=ceil(T/dt);
ppfac=-1.0; % pusher = 1, puller = -1.0
q0=1.0; % stresslet strength
U=1.0; %self swimming velocity;
b0=-(4.*pi.*d.^3).^(-1);

nstore=100;
t=0.0;
kappa=0.e2;
tol=0.1;

x_store=zeros([size(x_fd) ceil(Nt/nstore)]);
beta_store=zeros([size(x_fd) ceil(Nt/nstore)]);
x_store(:,:,1)=x_fd;
beta_store(:,:,1)=beta_fd;
j=1;
for n=1:Nt
    %compute forces due to steric forces
    f_s=stericforce(x_fd, kappa, tol);
    % compute velocity and dbeta based off point forces 
    u_s = stokesletvelocity(x_fd,x_fd,f_s,d);
    u_betas = stokesletbeta(x_fd,beta_fd,x_fd,f_s,d);
    % compute velocity and dbeta based off force doublets
    u_fd = fdvelocity(x_fd,beta_fd,d,q0,b0);
    u_betafd = fdbeta(x_fd,beta_fd,d,q0,b0);
    
    x_fd=x_fd+dt.*(u_s+u_fd);
    beta_fd=beta_fd+dt.*(u_betas+u_betafd);
    beta_fd=beta_fd./(vecnorm(beta_fd,2,2));
    

    if mod(n,nstore)==0
        j=j+1;
        x_store(:,:,j)=x_fd;
        beta_store(:,:,j)=beta_fd;
%         fname = sprintf('particles_%05d.vtk',j);
%         vtkwrite(fname, ...
%          'polydata', ...
%          'points', x_fd(:,1), x_fd(:,2), x_fd(:,3), ...
%          'vectors', 'orientation', beta_fd(:,1), beta_fd(:,2), beta_fd(:,3), ...
%          'vectors', 'force',       f_s(:,1), f_s(:,2), f_s(:,3), ...
%          'vectors', 'velocity_stokeslet',  u_s(:,1), u_s(:,2), u_s(:,3), ...
%          'vectors', 'velocity_fd',    u_fd(:,1), u_fd(:,2), u_fd(:,3));
    end
end

%%
if isfile('particles.h5')
    delete('particles.h5')
end

for k=1:nstore:Nt
    kk=round((k-1)/nstore)+1;
    coords = sprintf('/coords/step_%04d',kk);
    vel = sprintf('/orientation/step_%04d',kk);
    h5create('particles.h5',coords,[Nfd,3])
    h5write('particles.h5',coords,squeeze(x_store(:,:,kk)))
    h5create('particles.h5',vel,[Nfd,3])
    h5write('particles.h5',vel,squeeze(beta_store(:,:,kk)))
end

fileID = fopen('particles.xmf','w');
fprintf(fileID,'<?xml version="1.0" ?>\n');
fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fileID,'<Xdmf Version="3.0">\n');
fprintf(fileID,'<Domain>\n');
fprintf(fileID,'    <Grid Name="points" GridType="Collection" CollectionType="Temporal">\n');

for k=1:nstore:Nt
    kk=round((k-1)/nstore)+1;
    coords = sprintf('/coords/step_%04d',kk);
    vel = sprintf('/orientation/step_%04d',kk);
    fprintf(fileID,'        <Grid Name="step_%04d">\n',kk);
    fprintf(fileID,'            <Time Value="%.6f"/>\n', k*dt);
    fprintf(fileID,'            <Topology TopologyType="Polyvertex" NumberOfElements="%d">\n',Nfd);
    fprintf(fileID,'            </Topology>\n');
    fprintf(fileID,'            <Geometry GeometryType="XYZ">\n');
    fprintf(fileID,'                <DataItem Format="HDF" DataType="Float" Dimensions="3 %d">\n',Nfd);
    fprintf(fileID,'                    particles.h5:%s\n', coords);
    fprintf(fileID,'                </DataItem>\n');
    fprintf(fileID,'            </Geometry>\n');
    fprintf(fileID,'            <Attribute Name="Orientation" AttributeType="Vector" Center="Node">\n');
    fprintf(fileID,'                <DataItem Format="HDF" Dimensions="3 %d">\n',Nfd);
    fprintf(fileID,'                    particles.h5:%s\n',vel);
    fprintf(fileID,'                </DataItem>\n');
    fprintf(fileID,'            </Attribute>\n');
    fprintf(fileID,'        </Grid>\n');
end
fprintf(fileID,'    </Grid>\n');
fprintf(fileID,'</Domain>\n');
fprintf(fileID,'</Xdmf>\n');
fclose(fileID);

%% Plotting

for k=1:nstore:Nt
    kk=round((k-1)/nstore)+1;
    X = x_store(:,1,kk);
    Y = x_store(:,2,kk);
    Z = x_store(:,3,kk);
    
    U = beta_store(:,1,kk);
    V = beta_store(:,2,kk);
    W = beta_store(:,3,kk);
    
    figure(1)
    quiver3(X,Y,Z,U,V,W,'off','LineWidth',1.5)

    axis equal
    grid on

    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    xlim([-1,7])
    ylim([-3,3])
    zlim([-3,3])
    title(['Frame ',num2str(kk)])
    pause(0.1)
    drawnow


end


