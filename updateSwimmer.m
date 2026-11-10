function [u_s, beta_s, u_fd, beta_fd, u_r, beta_r] = updateSwimmer(i, x_list, beta_list, f_s, d, q0, b0, r0, Cr)
% UPDATESWIMMER calculates change in position and direction of force
% doublets for one agent
% Steric forces are assumed to already have been calculated, but this could
% be added here?
% Point sources of forces in f_s are assumed to be the other agents, whose
% positions are stored in x_fd
% - x_fd is a matrix of size Nx3, each row being an agent and each column
% one of the dimensions of the space
% - x_i is a row vector of length 3, giving the position of the agent of
% interest

% agent i
x_i = x_list(i,:);
beta_i = beta_list(i,:);

% measures of distance between agents
x = x_i - x_list;         % distance from agent i to all other agents
r = vecnorm(x,2,2);     % norm of distance
R = sqrt(d^2 + r.^2);   % regularized norm of distance

% dot products
dot_bx = x*beta_i';
dot_bf = f_s*beta_i';
dot_fx = dot(f_s,x,2);
dot_bb = beta_list*beta_i';
dot_xb = dot(beta_list, x, 2);

% vector triple products
cross_bxb = norm(beta_i)*x - dot_bx*beta_i;
cross_bfb = norm(beta_i)*f_s - dot_bf*beta_i;
cross_bbb = norm(beta_i)*beta_list - (beta_list*beta_i')*beta_i;
cross_bbx = dot_bx.*beta_list - dot_bb.*x;

% cross products
% cross_xb = cross(beta_list, x, 2);
cross_xb = beta_list(:,[2,3,1]).*x(:,[3,1,2]) - beta_list(:,[3,1,2]).*x(:,[2,3,1]);
cross_bbbx = cross_bbx(:,[2,3,1]).*beta_i([3,1,2]) - cross_bbx(:,[3,1,2]).*beta_i([2,3,1]);
cross_bb = beta_list(:,[2,3,1]).*beta_i([3,1,2]) - beta_list(:,[3,1,2]).*beta_i([2,3,1]);

% velocity due to steric forces
u_s = f_s.*(1 + d^2./R.^2) + dot(f_s,x,2).*x./R.^3;
u_s = sum(u_s./R,1)' ./ (8*pi);

% change in direction due to steric forces
beta_s = cross_bxb.*(dot_bf - 3*dot_bx.*dot_fx./R.^2) - cross_bfb.*dot_bx.*(1 + 3*d^2./R.^2);
beta_s = sum(beta_s./R.^3,1)'/(8*pi);

% velocity due to force doublets
u_fd = q0*(-x + 3*dot_bx.^2.*x./R.^2 + 3*dot_bx*beta_i*d^2./R.^2) + ...
    b0 * (-3*dot_bx.*x./R.^2 + (R.^2 - 3*d^2)*beta_i./R.^2);
u_fd = sum(u_fd./R.^3,1)'/(8*pi);

% change in direction due to force doublets
J1 = 6*dot_bb.*dot_xb + 3*(1 - 5*dot_xb.^2./R.^2).*dot_bx;
J2 = 3*d^2*(dot_bb - 5*dot_xb.*dot_bx./R.^2); % J2 and J3 use the same base vector, consider combining the two
J3 = -3*dot_bb + 15*dot_xb.*dot_bx./R.^2;
J4 = dot_bx.*(-3 + 15*d^2./R.^2);
beta_fd = (-q0*J1 + b0*J3).*cross_bxb + (-q0*J2 + b0*J4).*cross_bbb;
beta_fd = sum(beta_fd./R.^5,1)'/(8*pi);

% velocity due to rotlets
u_r = dot_xb.*cross_xb.*(6 + 15*d^2./R.^2) - Cr*cross_xb*15*d^3./R.^2;
u_r = r0*sum(u_r./R.^5,1)'/(8*pi);

% change in direction due to rotlets
J1 = cross_bbbx.*( (6 + 15*d^2./R.^2).*dot_bb - dot_bx.*(dot_xb.*(30./R.^2 + 105*d^2./R.^4) - Cr*105*d^3./R.^4));
J2 = cross_bb.*(dot_xb.*(6 + 15*d^2./R.^2) - Cr*15*d^3./R.^2);
beta_r = r0*sum((J1+J2)./R.^5,1)'/(8*pi);