function [t,z] = GlyT2_halt_ode(z0, k, kinv, c, tau, tspan, fix_GlyE, opts)

% scale
k = k  / c / tau;
kinv = kinv / tau;

Gly_e0 = z0(16);

cell_vol = 4/3 * pi * (1e-3 /2)^3;
bath_vol = 3e-4;
VR = cell_vol/bath_vol; % volume ratio of cell to surroundings
num_transporters = 1e11;


%Integration step
[t,z] = ode15s(@(t,z) odefcn(t,z, k, kinv, VR, num_transporters, fix_GlyE, Gly_e0), tspan, z0,opts);


%% set up ODE

function dzdt = odefcn(t,z, k, kinv, VR, num_transporters, fix_GlyE, Gly_e0)
dzdt = zeros(14,1);

y1 = z(1);
y2 = z(2);
y3 = z(3);
y4 = z(4);
y5 = z(5);

x1 =z(6);
x2 = z(7);
x3 = z(8);
x4 = z(9);
x5 = 1 - sum(z(1:9)); % conservation law

Na_i = z(11);
Na_e = z(12);

Cl_i =  z(13);
Cl_e = z(14);
Gly_i = z(15);
Gly_e = heaviside( tspan(end) / 2 - t) * z(16) +  heaviside(  t - tspan(end) / 2 ) * Gly_e0; 

Chi_i = Cl_i * Gly_i;
Chi_e = Cl_e * Gly_e;



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
k9 = k(9);
kinv9 = kinv(9);
k10 = k(10);
kinv10 = kinv(10);

dzdt(1) = k10 * x1  -(k1 *Chi_e +kinv10)*y1 + kinv1 *y2;
dzdt(2) =  k1 * Chi_e * y1  -(k2 * Na_e + kinv1 )  * y2 + kinv2 * y3;
dzdt(3) = k2* Na_e * y2 - (k3 * Na_e + kinv2)  * y3 +  kinv3 *  y4;
dzdt(4) = k3 * Na_e * y3 - (kinv3 + k4 * Na_e) * y4 + kinv4 * y5;
dzdt(5) = k4 * Na_e *y4 - (kinv4 + k5) * y5 +  kinv5 * x5;

dzdt(6) =  kinv10*y1  - (k10 + kinv9 * Chi_i) * x1 + k9 * x2;
dzdt(7) = kinv9*Chi_i * x1 -(k9 + kinv8 *Na_i )* x2 + k8 * x3;
dzdt(8) =  kinv8*Na_i * x2 -(k8  + kinv7 * Na_i )* x3 +  k7 * x4;
dzdt(9) =  kinv7 * Na_i * x3  - (k7 + kinv6 *Na_i) * x4 + k6 * x5;
dzdt(10) =   -sum(dzdt(1:9));



J_Na_i = (k6 *x5 - kinv6 * Na_i  * x4) + (k7 * x4 - kinv7 * Na_i  * x3) + (k8 * x3 - kinv8 * Na_i  * x2);
J_Na_e = (kinv2 * y3 - k2 *Na_e *y2) + (kinv3 *y4 - k3* Na_e * y3) + (kinv4 *y5 - k4* Na_e * y4);

% J_Cl = J_Gly
J_Chi_i=  k9 * x2 - kinv9 * Chi_i * x1;
J_Chi_e =  (kinv1 * y2 - k1* Chi_e* y1 );

dzdt(11) = J_Na_i; % intracellular sodium
dzdt(12) = 0;%J_Na_e; %extracellular sodium --fixed
dzdt(13) =  J_Chi_i; % intracellular chloride
dzdt(14) = 0; % extracellular chloride --fixed
dzdt(15) =  J_Chi_i; % intracellular Glycine
dzdt(16) = heaviside(  tspan(end) / 2 - t) * J_Chi_e; % extracellular Glycine
end
end