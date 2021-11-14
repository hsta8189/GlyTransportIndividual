function [k6, k7,k8,k9, kinv6, kinv7, kinv8, kinv9] =  GlyT2_bottom(n,c,g)

% define Rs here, node numbering is done in GLYT2_top.png
num_nodes = 21;
coop = 8;
R17 = 10^3 ;
R17inv = 300*coop;

R32inv = 100*coop;
R32 = 10^6;

R28 =10^6;
R28inv = 10^4 * coop;

R26 =10^4;
R26inv = 100*coop;

% node 2
R29 =10^3;
R29inv = 300*coop;

R30 =10^6;
R30inv = 10^4 * coop;

R31 = 10^4;
R31inv = 100*coop;

% node 3
R18 =10^6;
R18inv = 100*coop;

R20 = 10^6;
R20inv = 10^4 * coop;

R21 = 10^4;
R21inv = 100*coop;

%node 4
R19 = 10^3;
R19inv = 300*coop;

R27 = 10^6;
R27inv = 100*coop;

R29 = 10^3;
R29inv = 300*coop;

% node 5
R22 = 10^3;
R22inv = 300*coop;

R23 = 10^6;
R23inv = 10^4 * coop;


R25 = 10^6;
R25inv = 100*coop;

% node 6
R5 =10^6;
R5inv = 10^4*coop;

R13 =10^4;
R13inv = 100*coop;

% node 7
R6 =10^6;
R6inv = 100*coop;

R8  = 10^4;
R8inv = 100* coop;

% node 8
R7 = 10^3;
R7inv = 300*coop;

R14 = 10^4;
R14inv = 100*coop;

% node 9
R9 = 10^3;
R9inv = 300*coop;

R15 =10^6;
R15inv = 100*coop;

%node 10
R12 =10^3;
R12inv = 300*coop;

R16 =10^6;
R16inv = 10^4*coop;

% node 11;
R75 =10^5;
R75inv = 10^3 * coop;

R1 =10^6;
R1inv = 10^4 *coop;

% node 12
R3 = 10^6;
R3inv =100*coop;

R77 = 10^5;
R77inv = 10^3 * coop;

R10inv = 10^4 * coop;
R10 =10^6;

% node 13
R11 = 10^6;
R11inv = 100*coop;

% node 14
R4 = 10^6;
R4inv = 10^4 * coop;

R19 = 10^3;
R19inv = 300*coop;

% node 15
R1 = 10^6;
R1inv = 10^4 *coop;

R2 = 10^3;
R2inv = 300*coop;

% node 16
R67 = 10^4;
R67inv = 100;

%node 17
R68 = 10^3;
R68inv = 300;

%node 18
R69 = 10^5;
R69inv = 1000;

% node 19
R66 = 10^6;
R66inv = 10^4;

%node 20
R65 = 10^6;
R65inv =100;

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

rates(4,9) = R29 * n;
rates(9,4) = R29inv;

% node 5
rates(5,13) = R22 * n;
rates(13,5) = R22inv;

rates(5,9) = R23 * c;
rates(9,5) = R23inv;

rates(5,10) = R25 * g;
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

rates(8,15)  = R14 * n;
rates(15,8) = R14inv;

% node 9
rates(9,12) = R9 * n;
rates(12,9) = R9inv;

rates(9,15) = R15 * g;
rates(15,9) = R15inv;

%node 10
rates(10,14) = R12 * n;
rates(14, 10) = R12inv;

rates(10,15) = R16 * c;
rates(15,10) = R16inv;

% node 11;
rates(11,16) = R75 *n;
rates(16,11) = R75inv;

rates(11,18) = R1 * n;
rates(18,11) = R1inv;

% node 12
rates(12,18) = R3 * g;
rates(18,12) = R3inv;

rates(12,20) = R77 * n;
rates(20,12) = R77inv;

rates(12,13) = R10inv;
rates(13,12) = R10 * c;

% node 13
rates(13,14) = R11*g;
rates(14,13) = R11inv;

% node 14
rates(14,18) = R4 *c;
rates(18,14) = R4inv;

rates(14,19) = R19 * n;
rates(19,14) = R19inv;

% node 15
rates(15,17) = R1 * n;
rates(17,15) = R1inv;

rates(15, 18) = R2 * n;
rates(18,15) = R2inv;

% node 16
rates(16, 21) = R67 * n;
rates(21, 16) = R67inv;

%node 17
rates(17,21) = R68 * n;
rates(21,17) = R68inv;

%node 18
rates(18,21) = R69 * n;
rates(21,18) = R69inv;

% node 19
rates(19,21) = R66 * c;
rates(21,19) = R66inv;

%node 20
rates(20,21) = R65 * g;
rates(21,20) = R65inv;

%% Extra pathway for pseudo-equilibrium calculation
rates(2,21) = 300;
rates(21,2) = 300;

%% get k1 and k1inv
rates_2_8 = get_constants_ncoeffs(rates, [2, 8]);

k9 = rates_2_8(8,2);
kinv9 = rates_2_8(2,8);

rates_8_15 = get_constants_ncoeffs(rates, [8, 15]);
k8 = rates_8_15(15,8);
kinv8 = rates_8_15(8,15);

rates_15_16 = get_constants_ncoeffs(rates, [15, 16]);
k7 = rates_15_16(16,15);
kinv7 = rates_15_16(15,16);

rates_16_21 = get_constants_ncoeffs(rates, [16,  num_nodes]);
k6 = rates_16_21(21,16);
kinv6 = rates_16_21(16,21);



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
