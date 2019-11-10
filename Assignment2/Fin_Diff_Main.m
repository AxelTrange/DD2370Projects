%capacitor main -- Finite Difference main
%by Axel Trange

%Instructions:
%Objective: Reformulate the equations for Poisson Equation as a linear system
%           in matrix form. Solve it using x = A\B Matlab command.
%1) Check the slides, the Matlab code and Section 3.1.3 of the Textbook where 
%   we solved the problem using the Gauss-Seidel iteration that doesn't require
%   matrix formulation.
%2) Formulate the problem as matrix problem (it might be helpful to check 
%   the slides on finite differences).
%3) What is the equivalent of symmetry BC in matrix form?
%4) Plot the Potential on the grid and compare it with the result of the 
%   last iteration from the code we used in class.
%5) Perform a convergence test and plot it.
%6) Determine experimentally the order of truncation error

clear all; close all;

%Arguments:
a   = 1; %width of inner conductor
b   = 1; %height of inner conductor
c   =  2;%width of outer conductor
d   = 2; %height of outer conductor
n   = 20; %number of points in the x-direction (horizontal)
tol = 1E-3; %relative tolerance for capacitance
rel = 1; %relaxation parameter 

%Result of capacitance obtained in class:
cap_class = capacitor(a, b, c, d, n, tol, rel)

%Own result (matrix representation):
%cap_mat = capacitor_matrix(a, b, c, d, n, tol, rel)
%cap_mat_accurate = capacitor_matrix(a, b, c, d, 91, tol, rel) %uneven
cap_mat_accurate = capacitor_matrix(a, b, c, d, 90, tol, rel) %even

%Convergence test:
%nn = [3 5 15 25 35 45 55 65]; %NOT aligned for un-even
nn = [2 4 8 16 24 32 40 48 56 64]; %Aligned for even (due to rounding in capacitor.m)
%nn = [20];
cap_mat_list = zeros(1,numel(nn));
idx = 0;
for i = nn
    idx = idx+1;
    cap_mat = capacitor_matrix(a, b, c, d, i, tol, rel)
    cap_mat_list(idx) = cap_mat;
end
error_mat_list = -cap_mat_accurate + cap_mat_list;
refx = nn(1)*[1 1E2];
ref1 = error_mat_list(1)*[1 1E-2];
ref1_5 = error_mat_list(1)*[1 1E-3];
ref2 = error_mat_list(1)*[1 1E-4];
figure(3); loglog(nn,error_mat_list); grid on;
hold on; loglog(refx,ref1); loglog(refx,ref2); loglog(refx,ref1_5);
legend('Error','1st Order','2nd Order','3/2 Order')
xlim([2.5 100]);
xlabel('n'); ylabel('Error');


