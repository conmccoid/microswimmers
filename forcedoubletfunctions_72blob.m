function [H1,H2,J0,J1,J2,J3,J0por,J1por,J2por,J3por,H10]=forcedoubletfunctions_72blob(r,d)
% these correspond to the 7/2 blob.
% the stagnation point is located at:  s0 = - delta*p0/q0
%
% d = blob size
% r = distance from shilpalet location to the evaluation point (can be scalar or vector)
%
% Equations from Hoover, Boindala, and Cortez, 2025
R=(d.^2+r.^2).^(1/2);
H10=(1/4).*d.^(-1).*pi.^(-1);
H1=(1/4).*d.^2.*pi.^(-1).*R.^(-3)+(1/8).*pi.^(-1).*r.^2.*R.^(-3);
H2=(1/8).*pi.^(-1).*R.^(-3); % Eq 2.8, first line first term
J0=(-3/8).*pi.^(-1).*R.^(-5);% Eq 2.8, first line middle term
J1=(-3/8).*d.^2.*pi.^(-1).*R.^(-5);% Eq 2.8, first line last term
J2=(1/8).*pi.^(-1).*R.^(-5).*((-3).*d.^4+d.^2.*R.^2); %Eq 2.9, second term
J3=(-3/8).*d.^2.*pi.^(-1).*R.^(-5); % Eq 2.9, first term

J0por=(15/8).*pi.^(-1).*R.^(-7);
J1por=(15/8).*d.^2.*pi.^(-1).*R.^(-7);
J2por=(3/8).*pi.^(-1).*R.^(-7).*(5.*d.^4+(-1).*d.^2.*R.^2);
J3por=(15/8).*d.^2.*pi.^(-1).*R.^(-7);
end
