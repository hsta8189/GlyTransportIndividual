close all;
clearvars;

addpath('Calculate_constants')
addpath('Experimental')

% initial conditions
Na_i0 = 6e-3;
%Na_e0 = 150e-3;
q = 1.602e-19;
NC = 1e11; % num transporters in erdem et al
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
tmaxs = [1,1, 1, 1, 1];
fix_GlyE = false;
Na_exts = [ 150, 96, 10, 1] * 1e-3;
Na_ins = [6, 10, 10,10] * 1e-3;
qs = zeros(length(Na_exts), 1);
q_errs = zeros(length(Na_exts), 2);
linewidth=1.5;

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

debug = true; % whether to plot extra figures
legstring = [];
opts = odeset('RelTol',1e-10,'AbsTol',1e-14);
for i = 1:length(Na_exts)
    Na_e0 = Na_exts(i);
    Na_i0 = Na_ins(i);
    tspan = [0, tmaxs(i)];
    z0 = [y1, y2, y3, y4,y5, x1, x2, x3, x4,x5, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
    [~,z] = make_equilib_glyt2(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:10) = z(end, 1:10)';
    [t,z] = GlyT2_halt_ode(z0, k, kinv, c, tau, tspan, fix_GlyE, opts);
    

    Na_e = z(:, 12);
    I =  abs(q * 2/3 * NC * (k2 *Na_e .* z(:,2) - kinv2 *z(:,3) + k3 *Na_e .*  z(:,3) - kinv3 * z(:,4) + k4 *Na_e .*  z(:,4) - kinv4 * z(:,5)));
    
    figure(1);
    plot( t ,  z(:,16) * 1e6, 'LineWidth', linewidth )
    

        %[qs(i), q_errs(i, :)]  = get_hill_param(t,I);
         
         if i < 3
             figure(2);             
             plot( t(t < 0.5) ,  I(t < 0.5) / max(I(t < 0.5)),  'LineWidth', linewidth )
             figure(3); hold on;
             cutstep = 5e-6;
             plot( t(t > 0.5 + cutstep) ,  I(t > 0.5 +  cutstep) / max(I(t > 0.5 +  cutstep)),  'LineWidth', linewidth )
         end

    legstring = [ legstring sprintf("%d/%d", Na_e0 * 1e3, Na_i0 * 1e3)];
end
figure(2)
lgd = legend(legstring);
title(lgd, '[Na]_e/[Na]_i (mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.0es)',tau))
ylabel('I / I_{max}')
%xlim([0.018, max(t)])

figure(3)
lgd = legend(legstring);
title(lgd, '[Na]_e/[Na]_i (mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.0es)',tau))
ylabel('I / I_{max}')
xlim([0.5,0.5003])

figure(1)
lgd = legend(legstring);
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.2fs)',tau))
ylabel('[Gly]_e (\mu M)')

%exp_plot();
if ~debug
    % close figures used for debugging
    figure(4)
    close
    figure(5)
    close
end
function [q, err] = get_hill_param(ts, ys)
%


linewidth = 1.2;
if any(ys < 0)
    ys = ys / min(ys); % normalise
else
    ys = ys / max(ys); %normalise
end
if ys(end) - ys(end - 1) < 0 % need to invert the curve
    ys = max(ys) - ys; % shift to zero
else
    ys = ys - min(ys);
end
ys = ys / max(ys); % normalise again to account for all possible curve types

figure(5); hold on;
plot(ts, ys)

%snip off the data near the maximum that causes the hill curve to -> infinity.
% this makes the hill coefficient much more accurate (see figure 4)
ys = ys(ts < 0.7); 
ts = ts(ts < 0.7);

H = ys ./(1.1- ys);
ts = ts(15:end); % avoid log(0)
H = H(15:end);
figure(4); 
loglog(ts,H, 'LineWidth', linewidth)
hold on;
ylabel("I_{norm} / (1.1 - I_{norm})")
xlabel("Dimensionless Time")
f = fit(log(ts), log(H), 'poly1');
q = f.p1;
err = confint(f);
err = err(:,1);

end
