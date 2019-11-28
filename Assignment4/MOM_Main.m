
%MOM Assignment
%by Axel Trange and Itthinat Jongsuebchoke

clear all; close all;
format long

global r eta k0 gamma

%phi_inc = 0; %incoming angle
phi_inc = pi/4; %incoming angle

%r = 0.01; %1 cm radius
%r = 0.05; %5 cm radius
r = 0.20; %20 cm radius
f = 1E9; %frequency of inc. wave
c = 2.998E8; %speed of light
lambda = c/f; %wavelength
Approx_check = lambda/10; %Approximation
eta = 377; %free space impedance
k0 = 2*pi/lambda; %wave number
gamma = 1.781072418; %constant

Resol = [4 10 50 200 500 1000]; %Different resolutions
%Resol = [4]; %Different resolutions

Res_Accurate = 2000;
J_Accurate = MOM_Find_J(Res_Accurate,phi_inc); %Accurate solution
Len_surf = 2*pi*r;

%Error_list = zeros(1,numel(Resol));
itr = 0;
Strcat1 = string(missing);
for Res = Resol
    itr = itr+1;
    [J,x_n,y_n] = MOM_Find_J(Res,phi_inc);
    %J_dif = sum(abs(J_Accurate))/Res_Accurate - sum(abs(J))/Res;
    %Error_list(itr) = sqrt(sum(J_dif.*conj(J_dif))); %Pythagoras
    %Error_list(itr) = J_dif;
    %Error_list(itr) = J;
    J_mag = abs(J);
    J_phase = phase(J);

    figure(15);plot3(x_n,y_n,J_mag); grid on; 
    xlabel('x (m)');ylabel('y (m)');zlabel('magnitude'); hold on;
    figure(5);plot(linspace(0,360,numel(J_mag)),J_mag); grid on; 
    xlabel('Angle around contour (deg)');ylabel('magnitude'); hold on;
    Strcat1(itr) = strcat('n = ',num2str(Res)); %fixing str for legends
    
    figure(16);plot3(x_n,y_n,J_phase); grid on; 
    xlabel('x (m)');ylabel('y (m)');zlabel('phase (rad)'); hold on;
    figure(6);plot(linspace(0,360,numel(J_mag)),J_phase); grid on; 
    xlabel('Angle around contour (deg)');ylabel('phase (rad)'); hold on;
end
legend(Strcat1')
figure(5);legend(Strcat1')
figure(15);legend(Strcat1')
figure(16);legend(Strcat1')

%figure; loglog(Resol,abs(Error_list));xlabel('n');ylabel('Error');grid on;

phi = linspace(0,2*pi-(2*pi/Res),Res); %remove last point, duplicate
sigma_list = zeros(1,Res);
itr = 0;
for idx = phi
    itr = itr+1;
    J_Integral = sum(J.*exp(-1j*k0* (x_n'*cos(idx)+y_n'*sin(idx))))*2*pi*r/Res
    sigma_list(itr) = k0*eta^2/4 * abs(J_Integral)^2;
end

max(sigma_list) %display max RCS

figure;plot(phi*180/pi,sigma_list);xlabel('\Phi (deg)');ylabel('\sigma (m)');
hold on; grid on;

