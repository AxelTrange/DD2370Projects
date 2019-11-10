% --------------------------------------------------------------
% Compute capacitance per unit length of 
% a coaxial pair of rectangles
% --------------------------------------------------------------
% -- Provided code: Adjusted by Axel Trange
function cap = capacitor_matrix(a, b, c, d, n, tol, rel)

% Arguments:
%    a   =  width of inner conductor
%    b   = height of inner conductor
%    c   =  width of outer conductor
%    d   = height of outer conductor
%    n   = number of points in the x-direction (horizontal)
%    tol = relative tolerance for capacitance
%    rel = relaxation parameter 
%          (optimum is 2-c/n, where c is about pi)  
% Returns:
%    cap = capacitance per unit length [pF/m]

% Make grids
N = n+1;
h  = 0.5*c/n;                % Grid size
na = round(0.5*a/h);         % Number of segments on 'a'
x  = linspace(0,0.5*c,N);  % Grid points along x-axis
%n
m  = round(0.5*d/h);         % Number of segments on 'd'
M = m+1;
%h
mb = round(0.5*b/h);         % Number of segments on 'b'
y  = linspace(0,0.5*d,M);  % Grid points along y-axis


% --------------------------------------------------------------
%Matrix solution below:
A = diag(-4*ones(n^2,1)) + diag(ones(n^2-1,1),1) + diag(ones(n^2-1,1),-1); %contribution from x
A = A + diag(ones(n^2-n,1),n) + diag(ones(n^2-n,1),-n); %adding contribution from y
%C = zeros(n^2,1); %vector for f = A/c

%Fixing correction A_ij = 0 caused by boundary:
BC_skip_R = zeros(n^2-1,1);
BC_skip_L = zeros(n^2-1,1);
for nn = 1:n^2
    if mod(nn,n) == 0
        if nn ~= n^2
            BC_skip_R(nn) = -1;
            BC_skip_L(nn) = -1;
        end
    end
end
A = A + diag(BC_skip_R,1) + diag(BC_skip_L,-1); %correction

%Fixing boundary conditions:
%Dirichlet gives zero contribution (zero potential) - no change.
%Neumann gives contribution of same potential times nr of adjacent bounds.
%Inner conductor has potential zero - Add to source vector c.
for nn = 1:n^2
    if nn >= n^2 - (n-1) %bottom row, y
        A(nn,nn) = A(nn,nn)+1; %Neumann
    end
    if mod(nn,n) == 1 %left x
        A(nn,nn) = A(nn,nn)+1; %Neumann
    end
end

%Now; check what points are inside inner conductor for phi = 1:
C = zeros(n,m);  
for i = 1:na+1
  for j = 1:mb+1
    C(i,j)    = 1;
  end
end
C = flip(C); C_old = C;
C = C';
C = C(:);

%The only required equations inside the inner conductor is Phi = 1!
Idx = find(C); %non-zero indexes. Indexes of points inside inner conductor
Simple_row = zeros(1,n^2);
for i = Idx'
    Temp_row = Simple_row;
    Temp_row(i) = 1;
    A(i,:) = Temp_row;
end

%Potential answer:
f = A\C;
f = reshape(f,[n,n]);
f = f';

%Add a boundary in f corresponding to the zero-potential boundary:
%(This does not have to be part of the matrix)
New_f = zeros(n+1,n+1);
New_f((end-n+1):end,1:n) = f;
f = flip(New_f);

%Comment plot to run faster
figure(2); surf(f); xlabel('x');ylabel('y'); title('\Phi - Matrix');
cap = gauss(n,m,h,f);
% --------------------------------------------------------------


% --------------------------------------------------------------
% Compute capacitance from the potential
% --------------------------------------------------------------
function cap = gauss(n,m,h,f)

% Arguments:
%    n    = number of points in the x-direction (horizontal)
%    m    = number of points in the y-direction (vertical)
%    h    = cell size
%    f    = 2D-array with solution
% Returns:
%    cap = capacitance per unit length [pF/m]

q = 0;

for i = 1:n
  q = q + (f(i,m)+f(i+1,m))*0.5; %integrate along upper boundary
end

for j = 1:m
  q = q + (f(n,j)+f(n,j+1))*0.5; %integrate along right boundary
end

cap = q*4;           % 4 quadrants
cap = cap*8.854187;  % epsilon0*1e12 gives answer in pF/m
