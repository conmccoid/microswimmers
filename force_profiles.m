%% Force profiles
% This script measure the magnitude of forces as a function of
% (a) distance between agents
% (b) orientation between agents
% (c) number of agents (around a constant centre of mass)

clear
d = 1.0; % regularization parameter

ppfac = -1.0; % using puller
q0 = 4*pi*ppfac; % stresslet strength
U=1.0; % self swimming velocity (probably unnecessary for this test)
b0=-2*pi*d^3*U; % dipole (unclear?)
Cr=-1; % not sure what this is
r0=6*pi; % rotlet
kappa=1.e3; % steric force strength
tol=1*d; % tolerance for steric forces (likely to be the determining factor here)

%% Varying distance
dist_list = logspace(-1,4,100);
for ind = 1:length(dist_list)
    dist = dist_list(ind); % distance between test particles
    dir1 = pi/4;
    dir2 = 11*pi/6; % relative angle between test particles
    x_list = [0,0,0; dist, 0, 0];
    beta_list = [sin(dir1),cos(dir1),0;sin(dir2),cos(dir2),0];
    
    f_s = stericforce(x_list, kappa, tol);
    [u_s, beta_s, u_fd, beta_fd, u_r, beta_r] = updateSwimmer(1,x_list,beta_list,f_s, d, q0, b0, r0, Cr);
    mag_f(ind) = norm(f_s(1,:));
    mag_us(ind)= norm(u_s);
    mag_uf(ind)= norm(u_fd);
    mag_ur(ind)= norm(u_r);
    mag_bs(ind)= norm(beta_s);
    mag_bf(ind)= norm(beta_fd);
    mag_br(ind)= norm(beta_r);
end

%% Distance profiles
subplot(3,1,3)
semilogx(dist_list, mag_f,'k.--','linewidth',2)
xlabel('Distance between agents')
ylabel('Steric force magnitude')
subplot(3,1,2)
semilogx(dist_list, mag_us, 'r.--', dist_list, mag_uf, 'bo--', dist_list, mag_ur, 'm^--','linewidth',2)
ylabel('Change in velocity')
legend('Steric','Force doublets','Rotlets')
subplot(3,1,1)
semilogx(dist_list, mag_bs, 'r.--', dist_list, mag_bf, 'bo--', dist_list, mag_br, 'm^--','linewidth',2)
ylabel('Change in direction')
legend('Steric','Force doublets','Rotlets')

%% Varying direction
dir_list = linspace(0,2*pi, 100);
for ind1 = 1:length(dir_list)
    dist = 0.5; % distance between test particles
    dir1 = dir_list(ind1);
    x_list = [0,0,0; dist, 0, 0];
    beta_1 = [sin(dir1),cos(dir1),0];

    for ind2 = 1:length(dir_list)
        dir2 = dir_list(ind2);
        beta_2 = [sin(dir2),cos(dir2),0];
        
        f_s = stericforce(x_list, kappa, tol);
        [u_s, beta_s, u_fd, beta_fd, u_r, beta_r] = updateSwimmer(1,x_list,[beta_1;beta_2],f_s, d, q0, b0, r0, Cr);
        mag_f(ind1,ind2) = norm(f_s(1,:));
        mag_us(ind1,ind2)= norm(u_s);
        mag_uf(ind1,ind2)= norm(u_fd);
        mag_ur(ind1,ind2)= norm(u_r);
        mag_bs(ind1,ind2)= norm(beta_s);
        mag_bf(ind1,ind2)= norm(beta_fd);
        mag_br(ind1,ind2)= norm(beta_r);
    end
end

%% Direction profiles
% subplot(3,1,3)
% plot(dir_list, mag_f,'k.--','linewidth',2)
% xlabel('Angle between agents')
% ylabel('Steric force magnitude')
% subplot(3,1,2)
% plot(dir_list, mag_us, 'r.--', dir_list, mag_uf, 'bo--', dir_list, mag_ur, 'm^--','linewidth',2)
% ylabel('Change in velocity')
% legend('Steric','Force doublets','Rotlets')
% subplot(3,1,1)
% plot(dir_list, mag_bs, 'r.--', dir_list, mag_bf, 'bo--', dir_list, mag_br, 'm^--','linewidth',2)
% ylabel('Change in direction')
% legend('Steric','Force doublets','Rotlets')
subplot(2,3,2)
contourf(dir_list,dir_list,mag_uf)
xlabel('Angle 1')
ylabel('Angle 2')
title('u_{fd}')
subplot(2,3,3)
contourf(dir_list,dir_list,mag_ur)
title('u_r')
subplot(2,3,4)
contourf(dir_list,dir_list,mag_bs)
title('\beta_s')
subplot(2,3,5)
contourf(dir_list,dir_list,mag_bf)
title('\beta_{fd}')
subplot(2,3,6)
contourf(dir_list,dir_list,mag_br)
title('\beta_r')