function [k5, k6, k7,k5inv, k6inv, k7inv] =  GLYT1_bottom(n,c,g)
% define Rs here, node numbering is done in GLYT2_top.png
num_nodes = 16;
coop = 8;
R1 = 3*10^6;
R1inv = 10^3;

R2 = 3*10^6;
R2inv = 10^3;

R3 = 6*10^6;
R3inv = 200;

R4 = 10^6;
R4inv = 10^4;

R5 = 10^6;
R5inv = 10^4*coop;

R6 = 6*10^6;
R6inv = 200*coop;

R7 = 3*10^6;
R7inv = 10^3*coop;

R8 = 3*10^6;
R8inv = 10^3*coop;

R9 = 3*10^6;
R9inv = 10^3*coop;

R10 = 10^6;
R10inv = 10^4*coop;

R11 = 6*10^6;
R11inv = 200*coop;

R12 = 3*10^6;
R12inv = 10^3*coop;

R13 = 3*10^6;
R13inv = 10^3*coop;

R14 = 3*10^6;
R14inv = 10^3*coop;

R15 = 6*10^6;
R15inv = 200*coop;

R16 = 10^6;
R16inv = 10^4*coop;

R17 = 3*10^6;
R17inv = 10^3*coop;

R18 = 6*10^6;
R18inv = 200*coop;

R19 = 3*10^6;
R19inv = 10^3*coop;

R20 = 10^6;
R20inv = 10^4*coop;

R21 = 3*10^6;
R21inv = 10^3*coop;

R22 = 3*10^6;
R22inv = 10^3*coop;

R23 = 10^6;
R23inv = 10^4*coop;

R24 = 3*10^6;
R24inv = 10^3*coop;

R25 = 6*10^6;
R25inv = 200*coop;

R26 = 3*10^6;
R26inv = 10^3*coop;

R27 = 6*10^6;
R27inv = 200*coop;

R28 = 10^6;
R28inv = 10^4*coop;

R29 = 3*10^6;
R29inv = 10^3*coop;

R30 = 10^6;
R30inv = 10^4*coop;

R31 = 3*10^6;
R31inv = 10^3*coop;

R32 = 6*10^6;
R32inv = 200*coop;

%% set up matrix of rates
rates = zeros(num_nodes);

% I started at node one, and did forward + backwards rxn rates
% I then increased each node but never considered rates to lower node
% indexes to avoid repeating myself (not that it would matter)
% in the paper, the foward rate involves binding, so it is easy to keep
% track of which is forward and which is backward
% node 1
rates(1,6) = R17 * n;
rates(6,1) = R17inv;

rates(1,2) = R32inv;
rates(2,1) = R32 * g;

rates(1,8) = R28 * c;
rates(8,1) = R28inv;

rates(1,10) = R26 * n;
rates(10,1) = R26inv;

% node 2
rates(2,3) = R29 * n;
rates(3,2) = R29inv;

rates(2,4) = R30 *c;
rates(4,2) = R30inv;

rates(2,5) = R31 * n;
rates(5,2) = R31inv;

% node 3
rates(3,6) = R18 * n;
rates(6,3) = R18inv;

rates(3,7) = R20 * c;
rates(7,3) = R20inv;

rates(3,13) = R21 * n;
rates(13,3) = R21inv;

%node 4
rates(4,7) = R19 *n;
rates(7,4) = R19inv;

rates(4,8) = R27 * g;
rates(8,4) = R27inv;

rates(4,9) = R24 * n;
rates(9,4) = R24inv;

% node 5
rates(5,9) = R23 * c;
rates(9,5) = R23inv;

rates(5,13) = R22 * n;
rates(13,5) = R22inv;

rates(5,10) = R25 * n;
rates(10,5) = R25inv;

% node 6
rates(6,11) = R5 * c;
rates(11,6) = R5inv;

rates(6,14) = R13 *n;
rates(14,6) = R13inv;

% node 7
rates(7,11) = R6 * g;
rates(11,7) = R6inv;

rates(7,12) = R8 * n;
rates(12,7) = R8inv;

% node 8
rates(8,11) = R7 * n;
rates(11,8) = R7inv;

rates(8,15) = R14 * n;
rates(15,8) = R14inv;

% node 9
rates(9,12) = R9 * n;
rates(12,9) = R9inv;

rates(9,15) = R15 * g;
rates(15,9) = R15inv;

%node 10
rates(10,14) = R12 * n;
rates(14,10) = R12inv;

rates(10,15) = R16 * c;
rates(15,10) = R16inv;

% node 11;
rates(11,16) = R1 * n;
rates(16,11) = R1inv;

% node 12
rates(12,13) = R10inv;
rates(13,12) = R10 * c;

rates(12,16) = R3 * g;
rates(16,12) = R3inv;

% node 13
rates(13,14) = R11inv;
rates(14,13) = R11 * g;

% node 14
rates(14,16) = R4 *c;
rates(16,14) = R4inv;

% node 15
rates(15,16) = R2 * n;
rates(16,15) = R2inv;


rates(2,16) = 70;
rates(16,2) = 210;



%% get k1 and k1inv
rates_16_15 = get_constants_ncoeffs(rates, [16,15]);
k5 = rates_16_15(16,15);
k5inv = rates_16_15(15,16);

rates_8_15 = get_constants_ncoeffs(rates, [8,15]);
k6 = rates_8_15(15,8);
k6inv = rates_8_15(8,15);

rates_1_8 = get_constants_ncoeffs(rates, [1,8]);
k7 = rates_1_8(8,1);
k7inv = rates_1_8(1,8);
end


% names = string(1:21);
% zero_nodes = [];
% for i=1:21
% l(i) = nnz(rates2(i,:));
%     if l(i) == 0
%         zero_nodes = [zero_nodes i];
%     end
% end
% rates_db =double(subs(rates2, [n, c, g], [pi, pi, pi]));
% % names(zero_nodes) = [];
% % rates_db(zero_nodes,:) = [];
% 
% 
% G = graph(rates_db ~= 0);
% plot(G);
