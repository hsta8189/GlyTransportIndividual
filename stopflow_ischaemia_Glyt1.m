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
tmaxs = [1, 1, 1, 1, 1] ;
fix_GlyE = false;
Na_exts = [ 150, 96, 10, 1] * 1e-3;
Na_ins = [6, 10, 10,10] * 1e-3;
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

debug = true; % whether to plot extra figures
legstring = [];
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
for i = 1:length(Na_exts)
    Na_e0 = Na_exts(i);
    Na_i0 = Na_ins(i);
    tspan = [0, tmaxs(i)];
     z0 = [y1, y2, y3, y4, x1, x2, x3, x4, Na_i0, Na_e0, Cl_i0, Cl_e0, Gly_i0, Gly_e0]';
     [~,z] = make_equilib_glyt1(z0, k, kinv, c, tau, tspan / 10, true);
    
    z0(1:8) = z(end, 1:8);
    [t,z] = GlyT1_halt_ode(z0, k, kinv, c, tau, tspan, fix_GlyE, opts);
    
    z(t > tspan(end) / 2, 14) = Gly_e0;
    
    
    % plot slow changes on a different timescale
        
        timescale = t > 1.55e-5;
        
        Chi_i = z(:,11) .* z(:,13);
        I =  q  * NC * (0.7 * (k4 * z(:,4) - kinv4 *z(:,8)) + 0.3*(k7 * z(:,6) - kinv7 * Chi_i .*z(:,5)));
        


        % Due to accuracy of computer, once it reaches equilibrium it stops being a hill curve, cut off this section to get an accurate response
        % (perfect mathematical system => it wouldnt reach the eq. => hill
        % curve)

        %
        % the transporter config takes slightly longer to become resonable
        % when it is working badly

    
        figure(1); hold on;
         plot( t(timescale) ,  z(timescale,14) * 1e6,  'LineWidth', linewidth ) 
         figure(2); hold on;
        plot( t(timescale) , I(timescale) / max(abs(I(timescale))),  'LineWidth', linewidth ) 
    
    
    legstring = [ legstring sprintf("%d/%d", Na_e0 * 1e3, Na_i0 * 1e3)];
end
figure(1)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i(mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('[Gly]_e (\mu M)')

figure(2)
lgd = legend(legstring(:));
title(lgd, '[Na]_e/[Na]_i (mM)')
xlabel(sprintf('Dimensionless Time (t / \\tau), (\\tau = %.ds)',tau))
ylabel('I / I_{max}')
%ylim([0,5e-3])



if ~debug
    % close figures used for debugging
    figure(4)
    close
    figure(5)
    close
end

function [q, err] = get_hill_param(ts, ys,i)
%GET_HILL_PARAM Summary of this function goes here
%   Detailed explanation goes here


if any(ys < 0)
    ys = ys / min(ys); % normalise
else
    ys = ys / max(ys); %normalise
end
if ys(end) < ys(1) % need to invert the curve
    ys = 1 - ys; % shift to zero
else
    ys = ys - min(ys);
end

figure(5);
hold on;
ys = ys / max(ys);
plot(ts, ys);

%snip off the nonlinear tips (plot 5) that cause the hill coefficient to be
%underestimated
if i < 3
    ys = ys(ts < 0.12); 
    ts = ts(ts < 0.12);
elseif i == 4
    ys = ys(ts < 0.56); 
    ts = ts(ts < 0.56);
end


H = ys ./(1.1- ys);
ts = ts(15:end); % avoid log(0)
H = H(15:end);

figure(4); 
loglog(ts,H)
hold on;
ylabel("I_{norm} / (1.1 - I_{norm})")
xlabel("Dimensionless time")



f = fit(log(ts), log(H), 'poly1');
q = f.p1;
err = confint(f);
err = err(:,1);

end