% by Itthinat Jongsuebchoke
% This is a 1D FDTD simulation with pulse
% It displays a "movie" of the signal
% Size of the FDTD space
clear all;
close all;

ke=50;
% Position of the source
ks=ke/2;
% Number of time steps
nsteps=50;

% Cell size and time stepping
c0=3.e8;
dx=0.01;
% Use Courant condition
R = 0.5;
dt = R.*(dx/c0);
% Constants
cc=c0*dt/dx;
% Initialize vectors
E2z=zeros(nsteps,ke);   % E-field for 2nd order
Ez=zeros(1,ke);
Hy=zeros(1,ke);
% Gaussian pulse
t0=20;
spread=8;
%------ Start loop 1st order (Yee Scheme)
Ez=zeros(1,ke);
Hy=zeros(1,ke);
ExHist = [];
HyHist = [];
for t=1:nsteps
   % E field loop
   for k=2:ke-1
    Ez(k)=Ez(k)+cc*(Hy(k-1)-Hy(k));
   end
   % Source
   Ez(ks)=exp(-.5*((t-t0)/spread)^2);
   % H field loop
   for k=1:ke-1
      Hy(k)=Hy(k)+cc*(Ez(k)-Ez(k+1));
   end

   figure(1)
plot(Ez);axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');

%  pause()
   ExHist = [ExHist; Ez];
   HyHist = [HyHist; Hy];
end

%------ Start loop for 2nd order (Curl-Curl Eq.)
for t=3:nsteps
   % E field loop
   for k=2:ke-1
       E2z(t,k)=2*(1-cc^2)*E2z(t-1,k)-E2z(t-2,k)+cc^2*(E2z(t-1,k+1)+E2z(t-1,k-1));
   end
   % Source
   E2z(t,ks)=exp(-.5*((t-t0)/spread)^2);
end

figure(2);
plot(E2z(t,:));axis([1 ke -2 2]);
hold on
% plot(E1z(t,:))
legend('1st Order','2nd order');
xlabel('x');
ylabel('E_z');
hold off

%  pause()
end

%------ Dispersion ------
% spectrum1 = fft2(E1z);
spectrum1 = fft2(ExHist);
spectrum2 = fft2(E2z);

figure(10)
plot_spectrum1 = log(abs(spectrum1));
pcolor(plot_spectrum1(1:ke/2,1:nsteps/2))
shading interp

figure(11)
plot_spectrum2 = log(abs(spectrum2));
pcolor(plot_spectrum2(1:ke/2,1:nsteps/2))
shading interp

figure(12)
subplot(2,1,1)
plot_spectrum1 = log(abs(spectrum1));
pcolor(plot_spectrum1(1:ke/2,1:nsteps/2))
shading interp
title('1st Order')
xlabel('k\Deltax')
ylabel('\omega\Deltax/c')
subplot(2,1,2)
plot_spectrum2 = log(abs(spectrum2));
pcolor(plot_spectrum2(1:ke/2,1:nsteps/2))
shading interp
title('2nd Order')
xlabel('k\Deltax')
ylabel('\omega\Deltax/c')

%------ Comparison
figure(9)
subplot(2,2,1)
plot(E2z(20,:));axis([1 ke -2 2]);
hold on
plot(ExHist(20,:))
legend('1st Order','2nd order','south');
xlabel('x');
ylabel('E_z');
title('t=20');
grid on;
hold off
subplot(2,2,2)
plot(E2z(30,:));axis([1 ke -2 2]);
hold on
plot(ExHist(30,:))
legend('1st Order','2nd order','south');
xlabel('x');
ylabel('E_z');
title('t=30');
grid on;
hold off
subplot(2,2,3)
plot(E2z(40,:));axis([1 ke -2 2]);
hold on
plot(ExHist(40,:))
legend('1st Order','2nd order','south');
xlabel('x');
ylabel('E_z');
title('t=40');
grid on;
hold off
subplot(2,2,4)
plot(E2z(50,:));axis([1 ke -2 2]);
hold on
plot(ExHist(50,:))
legend('1st Order','2nd order','south');
xlabel('x');
ylabel('E_z');
title('t=50');
grid on;
hold off