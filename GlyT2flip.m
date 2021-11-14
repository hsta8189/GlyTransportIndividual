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
tau = 1e-3;
ntests = 100;
tmaxs = ones(1, ntests) ;
fix_GlyE = false;
Na_exts = linspace(0,150, ntests)* 1e-3;
Na_ins = linspace(0,150, ntests) * 1e-3;

Na_e_fixed = 150e-3;
Na_i_fixed = 10e-3;

qs = zeros(length(Na_exts), 1);
q_errs = zeros(length(Na_exts), 2);


load('data/GlyT2_ks.mat', 'k', 'kinv');
figure(1); hold on;
figure(2); hold on;

% get constants for calculating current
k2 = k(2);
kinv2 = kinv(2);
k3 = k(3);
kinv3 = kinv(3);
k4 = k(4);
kinv4 = kinv(4);
linewidth=1.5;

debug = false; % whether to plot extra figures
legstring = [];

Gly_de_finals = zeros(1, ntests);
Curr_de_finals = zeros(1, ntests);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
for i = 1:length(Na_exts)
    Na_e0 = Na_exts(i);
    Na_i0 = Na_i_fixed;
    tspan = [0, tmaxs(i)];
    
    % make better initial conditions
    z0 = [y1, y2, y3, y4,y5, x1, x2, x3, x4,x5, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
    [~,z] = make_equilib_glyt2(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:10) = z(end, 1:10);
    [t,z] = GlyT2_func(z0, k, kinv, c, tau, tspan, fix_GlyE, opts);
    
    
        
         Na_e = z(end, 12);
         Curr_de_finals(i) =  q * 2/3 * NC * (k2 *Na_e .* z(end,2) - kinv2 *z(end,3) + k3 *Na_e .*  z(end,3) - kinv3 * z(end,4) + k4 *Na_e .*  z(end,4) - kinv4 * z(end,5));
%         
        figure(1);
        plot( t(:),  z(:,16) * 1e6,  'LineWidth', linewidth ) 
        
        Gly_de_finals(i) = z(end, 16);
%         figure(2);
        
%         plot( t(timescale) ,  abs(I(timescale)) / max(abs(I)) ,  'LineWidth', linewidth ) 
        
    legstring = [ legstring sprintf("%d/%.0e", Na_e0 * 1e3, Na_i0 * 1e3)];
end
figure(1)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.0es)',tau))
ylabel('[Gly]_e (\mu M)')

figure(2)
plot(Na_exts * 1e3, Gly_de_finals * 1e6,  'LineWidth', linewidth );
xlabel('Extracellular Sodium (mM)')
ylabel('[Gly]_e final (\mu M)')

figure(3)
plot(Na_exts * 1e3, Curr_de_finals ,  'LineWidth', linewidth );
xlabel('Extracellular Sodium (mM)')
ylabel('I final (A)')




%% internal
Gly_di_finals = zeros(1, ntests);
Curr_di_finals = zeros(1, ntests);
legstring = [];

for i = 1:length(Na_exts)
    Na_e0 = Na_e_fixed;
    Na_i0 = Na_ins(i);
    tspan = [0, tmaxs(i)];
    
    
    % make better initial conditions
    z0 = [y1, y2, y3, y4,y5, x1, x2, x3, x4,x5, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
    [~,z] = make_equilib_glyt2(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:10) = z(end, 1:10);
    [t,z] = GlyT2_func(z0, k, kinv, c, tau, tspan, fix_GlyE, opts);
        
    
        
        Na_e = z(end, 12);
        Curr_di_finals(i) =  q * 2/3 * NC * (k2 *Na_e .* z(end,2) - kinv2 *z(end,3) + k3 *Na_e .*  z(end,3) - kinv3 * z(end,4) + k4 *Na_e .*  z(end,4) - kinv4 * z(end,5));

         
        figure(4); hold on;
        plot( t(:),  z(:,16) * 1e6,  'LineWidth', linewidth ) 
        
        Gly_di_finals(i) = z(end, 16);
%         figure(2);
        
%         plot( t(timescale) ,  abs(I(timescale)) / max(abs(I)) ,  'LineWidth', linewidth ) 
        
    legstring = [ legstring sprintf("%d/%.0e", Na_e0 * 1e3, Na_i0 * 1e3)];
end

figure(4)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('[Gly]_e (\mu M)')

figure(5)
plot(Na_ins * 1e3, Gly_di_finals * 1e6,  'LineWidth', linewidth );
xlabel('Intracellular Sodium (mM)')
ylabel('[Gly]_e final (\mu M)')

figure(6)
plot(Na_ins * 1e3, Curr_di_finals ,  'LineWidth', linewidth );
xlabel('Intracellular Sodium (mM)')
ylabel('I final (A)')




