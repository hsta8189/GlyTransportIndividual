function [t,z] = GlyT1_var_rates(z0, k, kinv, c, tau, tspan, fix_GlyE)

% scale
k = k  / c /  tau;
kinv = kinv / tau;

z0 = [z0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-14);

%Integration step
[t,z] = ode15s(@(t,z) odefcn(t,z), tspan, z0,opts);



%% set up ODE

function dzdt = odefcn(t,z)
dzdt = zeros(14,1);

y1 = z(1);
y2 = z(2);
y3 = z(3);
y4 = z(4);

x1 =z(5);
x2 = z(6);
x3 = z(7);
x4 = 1 - sum(z(1:7)); % conservation law

Na_i = z(9);
Na_e = z(10);

Cl_i =  z(11);
Cl_e = z(12);
Gly_i = z(13);
Gly_e = z(14);
Chi_i = Cl_i * Gly_i;
Chi_e = Cl_e * Gly_e;

if t == round(tspan(end)/2,0) % adjust constants near equilib
    path = 'temp/GlyT1.mat';
    make_GlyT1_consts(Na_i, Na_e, Cl_i, Cl_e, Gly_i, Gly_e, path)
    load(path, 'k', 'kinv')
end
    


k1 = k(1);
kinv1 = kinv(1);
k2 = k(2);
kinv2 = kinv(2);
k3 = k(3);
kinv3 = kinv(3);
k4 = k(4);
kinv4 = kinv(4);
k5 = k(5);
kinv5 = kinv(5);
k6 = k(6);
kinv6 = kinv(6);
k7 = k(7);
kinv7 = kinv(7);
k8 = k(8);
kinv8 = kinv(8);


dzdt(1) = k8 * x1  -(k1 *Chi_e +kinv8)*y1 + kinv1 *y2;
dzdt(2) =  k1 * Chi_e * y1  -(k2 * Na_e + kinv1 )  * y2 + kinv2 * y3;
dzdt(3) = k2* Na_e * y2 - (k3 * Na_e + kinv2)  * y3 +  kinv3 *  y4;
dzdt(4) = k3 * Na_e * y3 - (kinv3 + k4) * y4 + kinv4 * x4;

dzdt(5) = kinv8 * y1 - (k8 + kinv7*Chi_i) * x1 +  k7 * x2;
dzdt(6) = kinv7*Chi_i * x1 -(k7 + kinv6 *Na_i )* x2 + k6 * x3;
dzdt(7) =  kinv6*Na_i * x2 -(k6  + kinv5 * Na_i )* x3 +  k5 * x4;
dzdt(8) =  - sum(dzdt(1:7));


J_Na_i = (k5 *x4 - kinv5 * Na_i  * x3) + (k6 * x3 - kinv6 * Na_i  * x2);
J_Na_e = (kinv2 * y3 - k2 *Na_e *y2) - (kinv3 *y4 - k3* Na_e * y3);

% J_Cl = J_GLY
J_Chi_i=  k7 * x2 - kinv7 * Chi_i * x1;
J_Chi_e =  kinv1 * y2 - k1* Chi_e* y1;

dzdt(9) = J_Na_i; % intracellular sodium
dzdt(10) = 0;%J_Na_e; %extracellular sodium --fixed
dzdt(11) =  J_Chi_i; % intracellular chloride
dzdt(12) = 0; % extracellular chloride --fixed
dzdt(13) =  J_Chi_i; % intracellular glycine
if fix_GlyE % for flow type / replenising conditions
    dzdt(14) = 0;
else
    dzdt(14) =  J_Chi_e  ;  % extracellular Glycine
end
end
end