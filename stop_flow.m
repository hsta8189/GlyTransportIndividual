close all;
clearvars;

addpath('Calculate_constants')
addpath('Experimental')

% define initial conditions

Na_i0 = 6e-3;
%Na_e0 = 150e-3;
q = 1.602e-19;
NC = 1e11; % num transporters (for current calculation)
Cl_i0 = 9.4e-3;
Cl_e0 = 152e-3;

Gly_i0 = 2e-6;
Gly_e0 =  10e-6;


y1 = 0.5;
y2 = 0;
y3 = 0;
y4 = 0;
y5 = 0;

x1 = 0.5;
x2 = 0;
x3 = 0;
x4 = 0;
x5 = 0;


c = 1;
tau = 1e-1;
ntests = 5;
tmaxs = ones(1, ntests) ;
fix_GlyE = false;
Na_exts = linspace(0,150, ntests)* 1e-3;
Na_ins = linspace(0,150, ntests) * 1e-3;

Na_e_fixed = 50e-3;
Na_i_fixed = 10e-3;

qs = zeros(length(Na_exts), 1);
q_errs = zeros(length(Na_exts), 2);


load('data/GlyT1_ks.mat', 'k', 'kinv');
figure(1); hold on;
figure(2); hold on;

% get constants for calculating current
k7 = k(7);
kinv7 = kinv(7);
k4 = k(4);
kinv4 = kinv(4);
linewidth=1.5;

debug = false; % whether to plot extra figures
legstring = [];

Gly_de_finals = zeros(1, ntests);
Curr_de_finals = zeros(1, ntests);
J_de_finals = zeros(1, ntests);
for i = 1:length(Na_exts)
    Na_e0 = Na_exts(i);
    Na_i0 = Na_i_fixed;
    tspan = [0, tmaxs(i)];
    
    % make better initial conditions
    z0 = [y1, y2, y3, y4, x1, x2, x3, x4, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
    [~,z] = make_equilib_glyt1(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:8) = z(end, 1:8);
    [t,z] = GlyT1_halt_ode(z0, k, kinv, c, tau, tspan, fix_GlyE);
    
    
        
         Chi_i = z(:,11) .* z(:,13);
         Curr_de =  q  * NC * (0.7 * (k4 * z(:,4) - kinv4 *z(:,8)) + 0.3*(k7 * z(:,6) - kinv7 * Chi_i .*z(:,5)));
         
         Curr_de_finals(i) =  Curr_de(end);
         J_de_finals(i) =  q  * NC * (k4 * z(end,4) - kinv4 *z(end,8) );
%         
        figure(1);
        plot( t(:),  z(:,14) * 1e6,  'LineWidth', linewidth ) 
        
        figure(3); hold on;
        plot( t(:),  Curr_de,  'LineWidth', linewidth ) 
        
        figure(9); hold on;
        plot( t(:),  z(:,11) * 1e6,  'LineWidth', linewidth ) 
        
        Gly_de_finals(i) = z(end, 14);
%         figure(2);
        
%         plot( t(timescale) ,  abs(I(timescale)) / max(abs(I)) ,  'LineWidth', linewidth ) 
        
    legstring = [ legstring sprintf("%d/%.0e", Na_e0 * 1e3, Na_i0 * 1e3)];
end
figure(1)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('[Gly]_e (\mu M)')

figure(9)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('[Na]_i (\mu M)')

figure(2)
plot(Na_exts * 1e3, Gly_de_finals * 1e6,  'LineWidth', linewidth );
xlabel('Extracellular Sodium (mM)')
ylabel('[Gly]_e final (\mu M)')

figure(3);
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('I(A)')

figure(4); hold on;
plot(Na_exts * 1e3, Curr_de_finals ,  'LineWidth', linewidth );
xlabel('Extracellular Sodium (mM)')
ylabel('I final (A)')

plot(Na_exts * 1e3, J_de_finals ,  'LineWidth', linewidth );
legend(["Electrogenic Current", "Turnover Current"])








%% internal
Gly_di_finals = zeros(1, ntests);
Curr_di_finals = zeros(1, ntests);
J_di_finals = zeros(1, ntests);
legstring = [];
for i = 1:length(Na_exts)
    Na_e0 = Na_e_fixed;
    Na_i0 = Na_ins(i);
    tspan = [0, tmaxs(i)];
    
    
    % make better initial conditions
    z0 = [y1, y2, y3, y4, x1, x2, x3, x4, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
    [~,z] = make_equilib_glyt1(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:8) = z(end, 1:8);
    [t,z] = GlyT1_halt_ode(z0, k, kinv, c, tau, tspan, fix_GlyE);
        
    
        
         Chi_i = z(:,11) .* z(:,13);
         Curr_di =  q  * NC * (0.7 * (k4 * z(:,4) - kinv4 *z(:,8)) + 0.3*(k7 * z(:,6) - kinv7 * Chi_i .*z(:,5)));
         
         Curr_di_finals(i) =  Curr_di(end);
         J_di_finals(i) =  q  * NC * (k4 * z(end,4) - kinv4 *z(end,8) );
         
        figure(5); hold on;
        plot( t(:),  z(:,14) * 1e6,  'LineWidth', linewidth ) 
        
        figure(7); hold on;
        plot( t(:),  Curr_di ,  'LineWidth', linewidth ) 
        
        Gly_di_finals(i) = z(end, 14);
%         figure(2);
        
%         plot( t(timescale) ,  abs(I(timescale)) / max(abs(I)) ,  'LineWidth', linewidth ) 
        
    legstring = [ legstring sprintf("%d/%.0e", Na_e0 * 1e3, Na_i0 * 1e3)];
end

figure(5)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('[Gly]_e (\mu M)')



figure(6)
plot(Na_ins * 1e3, Gly_di_finals * 1e6,  'LineWidth', linewidth );
xlabel('Intracellular Sodium (mM)')
ylabel('[Gly]_e final (\mu M)')

figure(7); hold on;
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('I(A)')



figure(8); hold on;
plot(Na_ins * 1e3, Curr_di_finals ,  'LineWidth', linewidth );
xlabel('Intracellular Sodium (mM)')
ylabel('I final (A)')

plot(Na_ins * 1e3, J_di_finals ,  'LineWidth', linewidth );
legend(["Electrogenic Current", "Turnover Current"])






