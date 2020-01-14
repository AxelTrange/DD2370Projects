%Function to calculate matrices and get convergence:

function [A,B,no_ess,no_nat,A_nat,B_nat,Convergence, ...
    NrElements,Cutoff_FEM_20] = SolveMat_and_Conv(WG,PointsX,PointsY,Cutoff_Analyt_20)

global el2no no2xy noExt Lx Ly LxRidge LyRidge

%run mesh file to get el2no, no2xy, noInt, noExt
if strcmp(WG,'Rectangular')
    NrPtsX = PointsX;
    NrPtsY = PointsY;
    [el2no,no2xy,noInt,noExt] = Project_Mesh_Regular(NrPtsX,NrPtsY);
    %Cutoff_Analyt_20 = CutoffAnalytical()
elseif strcmp(WG,'Ridge')
    h_max = PointsX;
    [el2no,no2xy,noExt] = RidgeWaveguide(h_max,Lx,Ly,LxRidge,LyRidge)
    %Cutoff_Analyt_20 = RidgeWaveguide(h_max)
end

%Analytical cut-off solution:


% Initialization
%WG = 'rectangular';

mu0 = 4*pi*1e-7; 
c0 = 299792456;
eps0 = 1/(mu0*c0*c0);

% Reading the grid
%[no2xy, el2no, noExt] = ReadGrid([cross_section '_wg']);

noNum = size(no2xy,2);
elNum = size(el2no,2);

% Assembling
A = zeros(noNum);
B = zeros(noNum);

for elIdx = 1:elNum
  
  no = el2no(:,elIdx);  
  xy = no2xy(:,no);     
  
  [A_el, B_el] = CmpElMtx(xy);
  A(no,no) = A(no,no) + A_el;
  B(no,no) = B(no,no) + B_el;
  
end

% Getting the indices of the nodes.
no_ess = unique(noExt);
no_all = 1:noNum;
no_nat = setdiff(no_all, no_ess);

% Picking out the parts of the matrix and the vectors
% needed to solve the problem
A_ess    = A(no_nat,no_ess); %Corresponds to boundary
A_nat    = A(no_nat,no_nat); %Corresponds to inside

B_ess    = B(no_nat,no_ess);
B_nat    = B(no_nat,no_nat);



% Solving the eigenvalue problem
val_Hz = eig(A, B);
val_Ez = eig(A_nat, B_nat);

% Sorting and removing zero eigenvalue
val_Hz = sort(real(val_Hz));
val_Ez = sort(real(val_Ez));

val_Hz = val_Hz(2:length(val_Hz));

% Adding the eigenvalues together
val_all = [val_Hz; val_Ez];
tmte    = [zeros(size(val_Hz)); ones(size(val_Ez))];

[val_all, idx] = sort(val_all);
tmte = tmte(idx);

% Computing kt
kt_all = sqrt(val_all);

Cutoff_FEM_20 = c0*kt_all(1:20)/(2*pi);

% Plotting
%figure(26); clf;
figure; clf;
for idx = 1:20
  if (tmte(idx) == 0)
    % TE-mode
    plot(idx, c0*kt_all(idx)/(2*pi*1e9), 'o'), hold on
  else
    % TM-mode
    plot(idx, c0*kt_all(idx)/(2*pi*1e9), 's'), hold on
  end
end
xlabel('Mode number [-]')
ylabel('Cut-off frequency [GHz]')
title(['The ' WG ' waveguide (circles =' ...
       ' TE-modes and squares = TM-modes)']); grid on;
   
Diff = Cutoff_Analyt_20 - Cutoff_FEM_20

Convergence = abs(Diff);
%Convergence_list(:,PtIdx) = Convergence;
NrElements = elNum;
%NrElements_list(PtIdx) = NrElements;

end