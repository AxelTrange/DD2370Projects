%clear all
clear all; close all;

global el2no no2xy noExt Lx Ly LxRidge LyRidge

Lx = 2e-2;  % WG width
Ly = 1e-2;  % WG height
LxRidge = 1e-2; % Ridge width
LyRidge = 0.2e-2;  % Ridge heigth

%Points for rectangular mesh in X and Y. Lists need to be equal size!
%PointsY_List =  [3 4 5 7 10 12 15 17 20 24];
PointsY_List =  [3 5 9 15 18 24];
PointsX_List = 2*PointsY_List;

%Hmax = 1e-3*[5 4 3 2 1.5 1.3 1 0.75 0.6 0.5]; %Cell size (length) for ridge waveguide
Hmax = 1e-3*[1.5 1 0.75 0.5]; %Cell size (length) for ridge waveguide
Hmax_Analyt = 1E-3*0.45; %Running lower than 0.3 takes huge time.
%Hmax_Analyt = 1E-3*0.3; %Running lower than 0.3 takes huge time.

mu0 = 4*pi*1e-7; 
c0 = 299792456;
eps0 = 1/(mu0*c0*c0);

Convergence_Rect = [];
NrElements_Rect = [];
Cutoff_20_Rect_All = [];

Convergence_Ridge = [];
NrElements_Ridge = [];
Cutoff_20_Ridge_All = [];

%"Analytical" solutions:
Cutoff_Analyt_Rect = CutoffAnalytical();
[~,~,~,~,~,~,~,~,Cutoff_Analyt_Ridge] = SolveMat_and_Conv('Ridge',Hmax_Analyt,NaN,NaN);



WG = 'Ridge';

for PtIdx = 1:numel(Hmax)

    [A,B,no_ess,no_nat,A_nat,B_nat,Convergence, ...
        NrElements,Cutoff_20_Ridge] = SolveMat_and_Conv(WG,Hmax(PtIdx),NaN,Cutoff_Analyt_Ridge);
    
    Convergence_Ridge(:,PtIdx) = Convergence;
    NrElements_Ridge(PtIdx) = NrElements;
    Cutoff_20_Ridge_All(:,PtIdx) = Cutoff_20_Ridge;
    
end

WG = 'Rectangular';
for PtIdx = 1:numel(PointsX_List)

    [A,B,no_ess,no_nat,A_nat,B_nat,Convergence, ...
        NrElements,Cutoff_20_Rect] = SolveMat_and_Conv(WG,PointsX_List(PtIdx),...
        PointsY_List(PtIdx),Cutoff_Analyt_Rect);
    
    Convergence_Rect(:,PtIdx) = Convergence;
    NrElements_Rect(PtIdx) = NrElements;
    Cutoff_20_Rect_All(:,PtIdx) = Cutoff_20_Rect;
    
end

Conv_Mode1_Rect = Convergence_Rect(1,:);
Conv_Mode1_Ridge = Convergence_Ridge(1,:);

figure;loglog(NrElements_Rect,Conv_Mode1_Rect);grid on;
ylabel('\Delta f [Hz]'); xlabel('Nr of Elements');
Refx = NrElements_Rect(1) * [1 10^6];
Ref1y = Conv_Mode1_Rect(1) * [1 10^-6];
Ref2y = Conv_Mode1_Rect(1) * [1 10^(-2*6)];
Ref3y = Conv_Mode1_Rect(1) * [1 10^(-0.5*6)];
hold on; loglog(Refx,Ref1y); loglog(Refx,Ref2y); loglog(Refx,Ref3y);
xlim([NrElements_Rect(1) NrElements_Rect(end)]);
legend('Convergence 1st mode Rect','1st Order','2nd Order','1/2 Order')

figure;loglog(NrElements_Ridge,Conv_Mode1_Ridge);grid on;
ylabel('\Delta f [Hz]'); xlabel('Nr of Elements');
Refx = NrElements_Ridge(1) * [1 10^6];
Ref1y = Conv_Mode1_Ridge(1) * [1 10^-6];
Ref2y = Conv_Mode1_Ridge(1) * [1 10^(-2*6)];
Ref3y = Conv_Mode1_Ridge(1) * [1 10^(-0.5*6)];
hold on; loglog(Refx,Ref1y); loglog(Refx,Ref2y); loglog(Refx,Ref3y);
xlim([NrElements_Ridge(1) NrElements_Ridge(end)]);
legend('Convergence 1st mode Ridge','1st Order','2nd Order','1/2 Order')


%Extrapolation to zero-cell size for the cut-off freq of 1st mode (rect):
h_size_rect = Lx./PointsX_List;
h2 = h_size_rect.^2; %Convergence with respect to h is about ^2
pfit = polyfit(h2,Cutoff_20_Rect_All(1,:),3); %Choose m=3 here
pval = polyval(pfit,h2);
zero_size = polyval(pfit,0)
new_h = [h2 0];
new_pval = [pval zero_size];
figure;plot(new_h,new_pval,'r');grid on; %adding h=0 pot. to plot
hold on;plot(h2,Cutoff_20_Rect_All(1,:),'o');grid on;
ylabel('Cutoff frequency'); xlabel('h^2');
legend('3rd Order Extrapolation','Computed Frequency','Location','southeast')

%Extrapolation to zero-cell size for the cut-off freq of 1st mode (ridge):
h_size_ridge = Hmax;
h2r = h_size_ridge.^1.6 %about ^1.6 gives a linear plot in this case.
pfit = polyfit(h2r,Cutoff_20_Ridge_All(1,:),3); %Choose m=3 here
pval = polyval(pfit,h2r);
zero_size = polyval(pfit,0)
new_h = [h2r 0];
new_pval = [pval zero_size];
figure;plot(new_h,new_pval,'r');grid on; %adding h=0 pot. to plot
hold on;plot(h2r,Cutoff_20_Ridge_All(1,:),'o');grid on;
ylabel('Cutoff frequency'); xlabel('h^{1.6}');
legend('3rd Order Extrapolation','Computed Frequency','Location','southeast') 


%Rect to ridge cutoff frequency ratio:
f1f2_ratio_rect = Cutoff_20_Rect(2)/Cutoff_20_Rect(1)
f1f2_ratio_ridge = Cutoff_20_Ridge(2)/Cutoff_20_Ridge(1)

%the 20 first kt:
kt20_analyt_rect = Cutoff_Analyt_Rect/c0*2*pi;
kt20_rect = Cutoff_20_Rect/c0*2*pi;
kt20_ridge = Cutoff_20_Ridge/c0*2*pi;

%Visualize error plot for rectangular waveguide:
Error_Rect = Cutoff_Analyt_Rect - Cutoff_20_Rect;
figure;plot(abs(Error_Rect)); hold on; grid on;
xlabel('Mode number');ylabel('\Delta f [Hz]')


%For visualizing:
%From here it is only a matter of displaying what we got.
%We got the fields, the elements numbers, the coordinates.

%H here corresponds to Hz, giving Transversal E, i.e. TE Modes!
%First seen mode is TE10, then TE20, TE01, etc.
[Hz,ktHz] = eig(A,B);
ktHz = diag(ktHz);
[ktHz,Idx] = sort(real(ktHz));
ktHz = ktHz(2:length(ktHz));
Hz = Hz(:,Idx);
H_Full = Hz(:,1); %First is TE00, so Hz(:,1) is irrelevant!

%E here corresponds to Ez, giving Transversal H, i.e. TM Modes!
%First seen mode is TM11, then TM21, etc.!
[Ez,ktEz] = eig(A_nat,B_nat);
ktEz = diag(ktEz);
[ktEz,sort_Idx] = sort(real(ktEz));
ktEz = ktEz(2:length(ktEz));
Ez = Ez(:,sort_Idx);
fc_eig_20 = sqrt(ktEz(1:20))*c0/(2*pi);

%Visualize longitudinal fields:
E_Full = zeros(size(no_ess,1)+size(no_nat,2),1);
for idx = 1:5 %Visualize first 5 TM modes
    E_Full(no_nat) = Ez(:,idx);
    VisualizeMode(no2xy, el2no, E_Full,'E_z');
end

for idx = 1:5 %Visualize first 5 TE modes
    H_Full = Hz(:,idx);
    VisualizeMode(no2xy,el2no,H_Full,'H_z');
end


%Visualize curl of longitudinal fields (transversal fields):
E_FullX = zeros(size(no_ess,1)+size(no_nat,2),1);
for Mode = 1:5 %Visualize curl of 5 first TM modes
    E_FullX(no_nat) = Ez(:,Mode);
    Visualize_Curl(no2xy,el2no,Mode,E_FullX,'H_t (TM)')
end

for Mode = 1:5 %Visualize curl of 5 first TE modes
    Visualize_Curl(no2xy,el2no,Mode,Hz(:,Mode),'E_t (TE)')
end




