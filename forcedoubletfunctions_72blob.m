function [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d)
% these correspond to the 7/2 blob.
% the stagnation point is located at:  s0 = - delta*p0/q0
%
% d = blob size
% r = distance from shilpalet location to the evaluation point (can be scalar or vector)
%
R=(d.^2+r.^2).^(1/2);
H10=(1/4).*d.^(-1).*pi.^(-1);
H1=(1/4).*d.^2.*pi.^(-1).*(d.^2+r.^2).^(-3/2)+(1/8).*pi.^(-1).*r.^2.*(d.^2+r.^2).^(-3/2);
H2=(1/8).*pi.^(-1).*(d.^2+r.^2).^(-3/2);
J0=(-3/8).*pi.^(-1).*R.^(-5);
J1=(-3/8).*d.^2.*pi.^(-1).*R.^(-5);
J2=(1/8).*pi.^(-1).*R.^(-5).*((-3).*d.^4+d.^2.*R.^2);
J3=(-3/8).*d.^2.*pi.^(-1).*R.^(-5);
J0por=(15/8).*pi.^(-1).*R.^(-7);
J1por=(15/8).*d.^2.*pi.^(-1).*R.^(-7);
J2por=(3/8).*pi.^(-1).*R.^(-7).*(5.*d.^4+(-1).*d.^2.*R.^2);
J3por=(15/8).*d.^2.*pi.^(-1).*R.^(-7);

% R  =sqrt(d.^2 + r.^2);
% H10=1./(4*d*pi);
% H1 =(1+r.^2/2)./(4*d.^2*pi*R.^3);
% H2 =1./(8*pi*R.^3);
% J0 =-3./(8*pi*R.^5);
% J1 =-3*d.^2./(8*pi*R.^5);
% J2 =(-3*d.^4+d.^2.*R.^2)./(8*pi*R.^5);
% J3 =(-3*d.^2)./(8*pi*R.^5);
% J0por=15./(8*pi*R.^7);
% J1por=15
end
