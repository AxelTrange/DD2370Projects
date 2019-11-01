%by Axel Trange and Itthinat Jongsuebchoke
%Simpson rule convergence test: Assignment 1 Computational Electromag.
clear all;
close all;

rule = 'simpson'; %method
%rule = 'midpoint'; %method
z = 1; %height in m
a = 1; %square side in m
nn = [5 7 10 15 20]; %grid spacing h = a/n
%nn = [5 50 500 1200 2000]; %grid spacing h = a/n
%nn = [20 40 60 80 100]; %grid spacing h = a/n

pot_list = zeros(1,numel(nn));
time_list = zeros(1,numel(nn));
stats_list = zeros(1,numel(nn));

i = 1;
for n = nn
    
    pot = integr(z,a,n,rule);
    pot_list(i) = pot;
   
    stats = profile('INFO');
    struct_integr = stats.FunctionTable(2); %integr in channel 2.
                                            %see profile viewer.
    time_list(i) = struct_integr.TotalTime;
    
    profile viewer
   
    i = i+1;
end

Sol_Accurate = integr(z,a,2000,rule);
Error_list = abs(Sol_Accurate - pot_list);

ref_x = [1E1 1E3];
ref_o1_y = [1E-6 1E-8];
ref_o2_y = [1E-6 1E-10];
ref_o3_y = [1E-6 1E-12];
ref_o4_y = [1E-6 1E-14];

h = (a./nn); %grid spacing
if strcmp(rule,'midpoint')
    h2 = h.^2; %Power 2 for midpoint
end
if strcmp(rule,'simpson')
    h2 = h.^4; %Power 4 for simpson
end
%h2 = h.^4;
figure;plot(h2,pot_list,'o');grid on; ylabel('\Phi');
if strcmp(rule,'midpoint')
    xlabel('h^2');
end
if strcmp(rule,'simpson')
    xlabel('h^4');
end
figure;plot(nn,time_list,'-o');grid on; xlabel('n'); ylabel('time (s)');
figure;loglog(nn,Error_list,'o');grid on; xlabel('n'); ylabel('Error');
hold on; loglog(ref_x,ref_o1_y); loglog(ref_x,ref_o2_y);
loglog(ref_x,ref_o3_y); loglog(ref_x,ref_o4_y);
legend('Computed Error','Order 1','Order 2','Order 3', 'Order 4');
%-> Order of convergence seems to be Order 4 (Simpson)

pfit = polyfit(h2,pot_list,2); %Choose m=2 here
pval = polyval(pfit,h2);
pot_list
pval
zero_size = polyval(pfit,0)
new_h = [h2 0];
new_pval = [pval zero_size];
figure;plot(new_h,new_pval,'r');grid on; %adding h=0 pot. to plot
hold on;plot(h2,pot_list,'o');grid on; ylabel('\Phi');
if strcmp(rule,'midpoint')
    xlabel('h^2');
end
if strcmp(rule,'simpson')
    xlabel('h^4');
end
legend('2nd Order Extrapolation','Computed Potential','Location','southeast')

