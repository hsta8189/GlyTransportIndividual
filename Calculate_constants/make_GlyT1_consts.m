function make_GlyT1_consts(Na_i, Na_e, Cl_i, Cl_e, Gly_i, Gly_e, path)
    addpath("Calculate_constants")
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
    

    [k5, k6, k7, kinv5, kinv6, kinv7] =  GlyT1_bottom(Na_i, Cl_i, Gly_i);
    [k1, k2,k3, kinv1, kinv2, kinv3] =  GlyT1_top(Na_e, Cl_e,Gly_e);


    k4 = 70;
    kinv4 = 70;
    k8 = 210;
    kinv8 = 210;

    % voltage dependence of rate constants
    V = -60 *1e-3; % 60 mV, fixed;
    kB = physconst('Boltzmann');
    T = 300; % temp in kelvin
    q = 1.602e-19; % fundamental charge in coulombs

    K1 = kinv1/k1;
    K2 = kinv2 / k2;
    K3 = kinv3 / k3;
    K5 = kinv5/ k5;
    K6 = kinv6/k6;
    K8 = kinv8/ k8;


    kappa = 1 /(K1 *K2 * K3 *  K5 * K6  * K8 );

    % K1 is voltage dependent,  all Kis are
    % constant, so detailed balance requires:
    kinv7 = k7 * kappa *exp( 0.3 *q * V / (kB * T)); 
    kinv4 = k4 * kappa *exp( 0.7* q * V / ( kB * T)); 

    k = [k1 k2 k3 k4 k5 k6 k7 k8];
    kinv = [kinv1 kinv2 kinv3 kinv4 kinv5 kinv6 kinv7 kinv8];
    
    
    save(path, 'k', 'kinv');

end