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
tau = 1;
ntests = 20;
tmaxs = ones(1, ntests) ;
fix_GlyE = false;
Na_exts = linspace(0,150, ntests)* 1e-3;
Na_ins = zeros(1, ntests);
Gly_es = [1,10,100]*1e-6; %Gly_es = [1, 10,20,50,100,1000] * 1e-6;

Na_e_fixed = 50e-3;
Na_i_max = 300e-3;


load('data/GlyT1_ks.mat', 'k', 'kinv');
figure(1); hold on;

% get constants for calculating current
k7 = k(7);
kinv7 = kinv(7);
k4 = k(4);
kinv4 = kinv(4);
linewidth=1.5;

debug = false; % whether to plot extra figures
legstring = [];
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
for Gly_e0 = Gly_es
    for i = 1:length(Na_exts)
        Na_e0 = Na_exts(i);
        Na_i_inc = Na_e0/10;
        Na_i_inc = max(Na_i_inc, 0.1e-3);
        Na_i0 = Na_e0 - Na_i_inc;
        tspan = [0, tmaxs(i)];
        Gly_e_ss = 0;

        while Gly_e_ss < Gly_e0 && Na_i0 <= Na_i_max
             Na_i0 = Na_i0 + Na_i_inc;
            % make better initial conditions
            z0 = [y1, y2, y3, y4, x1, x2, x3, x4, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
            [~,z] = make_equilib_glyt1(z0, k, kinv, c, tau, tspan / 10, true);

            % calculate timeseries
            z0(1:8) = z(end, 1:8);
            [t,z] = GlyT1_func(z0, k, kinv, c, tau, tspan, fix_GlyE, opts);

            % get steady state
            Gly_e_ss = z(end,14);
        end
        Na_ins(i) = Na_i0;
    end
    figure(1);
    plot(Na_exts(Na_ins <= 295e-3) * 1e3, Na_ins(Na_ins <= 295e-3) * 1e3,  'LineWidth', linewidth );
    legstring = [legstring, sprintf("%d", round(Gly_e0*1e6,0))];
end
lgd = legend(legstring);
title(lgd, "Extracellular Glycine (\muM)")
xlabel('Extracellular Sodium (mM)')
ylabel('Intracellular Sodium (mM)')

f = fit(Na_exts', Na_ins', 'poly1');
