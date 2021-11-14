function make_GlyT2_consts(Na_i, Na_e, Cl_i, Cl_e, Gly_i, Gly_e)

    if nargin == 0
            % rates in erdem et all were measured under typical cell conditions, so it
        % makes the most sense to evalutate rate constants under those conditions

        % typical conditions from erdem et al - not initial conditions,
        % these are only for evaulating k1 and k1inv
        Na_i = 5.9e-3;
        Na_e = 150e-3;

        Cl_i = 9.3e-3;
        Cl_e = 154e-3;

        Gly_i =  2.2e-6;
        Gly_e = 1e-7;
    end
    

    [k6, k7,k8,k9, kinv6, kinv7, kinv8, kinv9] =  GlyT2_bottom(Na_i, Cl_i, Gly_i);
    [k1, k2,k3,k4, kinv1, kinv2, kinv3, kinv4] =  GlyT2_top(Na_e, Cl_e,Gly_e);


    k5 = 300;
    kinv5 = 300;
    k10 = 300;
    kinv10 = 300;

    % voltage dependence of rate constants
    V = -60 *1e-3; % 60 mV, fixed;
    kB = physconst('Boltzmann');
    T = 300; % temp in kelvin
    q = 1.602e-19; % fundamental charge in coulombs

    K1 = kinv1/k1;
    K5 = kinv5/ k5;
    K6 = kinv6/k6;
    K7 = kinv7/ k7;
    K8 = kinv8/ k8;
    K9 = kinv9/k9;
    K10 = kinv10/k10;

    kappa = 1 /(K1 * K5 * K6 * K7 * K8 * K9 * K10);

    % K1 is voltage dependent,  all Kis are
    % constant, so detailed balance requires:
    kinv2 = k2 * kappa *exp( 2* q * V / (3 * kB * T)); 
    kinv3 = k3 * kappa *exp( 2* q * V / (3 * kB * T)); 
    kinv4 = k4 * kappa *exp( 2* q * V / (3 * kB * T)); 

    k = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10];
    kinv = [kinv1 kinv2 kinv3 kinv4 kinv5 kinv6 kinv7 kinv8 kinv9 kinv10];
    
    save('GlyT2_ks.mat', 'k', 'kinv');

end