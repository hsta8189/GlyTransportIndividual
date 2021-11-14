close all;
clearvars;

addpath('Calculate_constants')
addpath('Experimental')
q = 1.602e-19;
NC = 1e11; % num transporters in erdem et al

Na_i0 = 10e-3;
Na_e0 = 150e-3;

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
tspan = [0,1];
fix_GlyE = false;
nGlys = 100;
Nas = [150,1] * 1e-3;
Gly_exts = linspace(0,300, nGlys) *1e-6  ;
load('data/GlyT2_ks.mat', 'k', 'kinv');
% get constants for calculating current
k7 = k(7);
kinv7 = kinv(7);
k4 = k(4);
kinv4 = kinv(4);

j = 1;
for Na_e0 = Nas
    figure(j);
    
    hold on;
    I_ss = zeros(1, nGlys); % steady state current;
    for i = 1:nGlys
        Gly_e0 = Gly_exts(i);
        z0 = [y1, y2, y3, y4, y5, x1, x2, x3, x4, x5, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
        [t,z] = GlyT2_func(z0, k, kinv, c, tau, tspan,fix_GlyE);

        Na_e = z(:,12);
        Chi_i = z(:,11) .* z(:,13);
        I =  abs(q  * NC * (0.7 * (k4 * z(:,4) - kinv4 *z(:,8)) + 0.3*(k7 * z(:,6) - kinv7 * Chi_i .*z(:,5))));
        figure(j)
        plot( t ,  I )
        I_ss(i) = I(end);

        %legstring = [ legstring sprintf("%d/%d", Na_e0 * 1e3, Na_i0 * 1e3)];
    end

    %title(lgd, '[Na]_e/[Na]_i (mM)')
    xlabel('dimensionless time')
    ylabel('I')

    figure(3); hold on;
    plot(Gly_exts * 1e6, 1- I_ss / max(I_ss) , 'LineWidth', 1.5);
    j= j+1;
end

lgd = legend(["150/10","1/10"]);
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel("Extracellular Glycine concentration (\mu M)")
ylabel("I_{steady state} / I_{steady state, max}")