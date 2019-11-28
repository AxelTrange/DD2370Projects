%Function to find the Z-matrix, plot it, and to also solve the
%Method of Moments problem in finding the current density J.

%by Axel Trange and Itthinat Jongsuebchoke

function [J,x_n,y_n] = MOM_Find_J(Resol, phi_inc)
%
global r eta k0 gamma
E0 = 1;
w_n = 2*pi*r/Resol;
phi = linspace(0,2*pi-(2*pi/Resol),Resol); %remove last point, duplicate

H_02 = @(x) besselh(0,2,x);

x_n = r*cos(phi);
y_n = r*sin(phi);
x_mat = x_n - x_n';
y_mat = y_n - y_n';

R_mn = sqrt((x_mat).^2 + (y_mat).^2);

Z_nn = k0*eta*w_n/4 * (1-1j*2/pi * (log(gamma*k0*w_n/4)-1));
Z_mn = k0*eta*w_n/4 * H_02(k0*R_mn);
b_m = E0*exp(-1j*k0*(x_n'*cos(phi_inc)+y_n'*sin(phi_inc)));

Z_mn(boolean(eye(Resol))) = Z_nn;
figure(10);imagesc(abs(Z_mn))

J = Z_mn\b_m;

end
