%close all
%clear all

function Cutoff_20 = CutoffAnalytical()

% Parameters
c0 = 299792456;    % Speed of light in free space
Lx = 2e-2;  % WG width
Ly = 1e-2;  % WG height
NrModes = 20;
 
nnx = 5;
nny = 5;

fc_listTE = [];
fc_listTM = [];

count=0;
for nx = 1:nnx
    for ny = 1:nny
       count=count+1;
       fc_listTE(count) = cutoff(nx-1,ny-1,Lx,Ly,c0);
       fc_listTM(count) = cutoff(nx,ny,Lx,Ly,c0);
    end
end

[out,idx] = sort([fc_listTE fc_listTM]);

Cutoff_20 = out(2:NrModes+1)';

figure(28);plot(Cutoff_20/1E9,'o');grid on;
xlabel('Mode number'); ylabel('Cut-off frequency [GHz]');

function fc = cutoff(nx,ny,Lx,Ly,c0)
fc = c0*(1/(2*pi))*sqrt((pi*nx/Lx)^2+(pi*ny/Ly)^2);
end

end
