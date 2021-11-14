function [k1, k2, k3,k1inv, k2inv, k3inv] =  GLYT1_top(n,c,g)

% define Rs here, node numbering is done in GLYT2_top.png
num_nodes = 16;
coop = 8;

R33 = 3*10^6;
R33inv = 10^3;

R34 = 3*10^6;
R34inv = 10^3;

R35 = 6*10^6;
R35inv = 200;

R36 = 10^6;
R36inv = 10^4;

R37 = 10^6;
R37inv = 10^4*coop;

R38 = 6*10^6;
R38inv = 200*coop;

R39 = 3*10^6;
R39inv = 10^3*coop;

R40 = 3*10^6;
R40inv = 10^3*coop;

R41 = 3*10^6;
R41inv = 10^3*coop;

R42 = 10^6;
R42inv = 10^4*coop;

R43 = 6*10^6;
R43inv = 200*coop;

R44 = 3*10^6;
R44inv = 10^3*coop;

R45 = 3*10^6;
R45inv = 10^3*coop;

R46 = 3*10^6;
R46inv = 10^3*coop;

R47 = 6*10^6;
R47inv = 200*coop;

R48 = 10^6;
R48inv = 10^4*coop;

R49 = 3*10^6;
R49inv = 10^3*coop;

R50 = 6*10^6;
R50inv = 200*coop;

R51 = 3*10^6;
R51inv = 10^3*coop;

R52 = 10^6;
R52inv = 10^4*coop;

R53 = 3*10^6;
R53inv = 10^3*coop;

R54 = 3*10^6;
R54inv = 10^3*coop;

R55 = 10^6;
R55inv = 10^4*coop;

R56 = 3*10^6;
R56inv = 10^3*coop;

R57 = 6*10^6;
R57inv = 200*coop;

R58 = 3*10^6;
R58inv = 10^3*coop;

R59 = 6*10^6;
R59inv = 200*coop;

R60 = 10^6;
R60inv = 10^4*coop;

R61 = 3*10^6;
R61inv = 10^3*coop;

R62 = 10^6;
R62inv = 10^4*coop;

R63 = 3*10^6;
R63inv = 10^3*coop;

R64 = 6*10^6;
R64inv = 200*coop;

R65 = 210;
R65inv = 210;

R66 = 70;
R66inv = 70;
%% set up matrix of rates
rates = zeros(num_nodes);

% I started at node one, and did forward + backwards rxn rates
% I then increased each node but never considered rates to lower node
% indexes to avoid repeating myself (not that it would matter)
% in the paper, the foward rate involves binding, so it is easy to keep
% track of which is forward and which is backward
% node 1
rates(1,6) = R49 * n;
rates(6,1) = R49inv;

rates(1,2) = R64inv;
rates(2,1) = R64 * g;

rates(1,8) = R60 * c;
rates(8,1) = R60inv;

rates(1,10) = R58 * n;
rates(10,1) = R58inv;

% node 2
rates(2,3) = R61 * n;
rates(3,2) = R61inv;

rates(2,4) = R62 *c;
rates(4,2) = R62inv;

rates(2,5) = R63 * n;
rates(5,2) = R63inv;

% node 3
rates(3,6) = R50 * n;
rates(6,3) = R50inv;

rates(3,7) = R52 * c;
rates(7,3) = R52inv;

rates(3,13) = R53 * n;
rates(13,3) = R53inv;

%node 4
rates(4,7) = R51 *n;
rates(7,4) = R51inv;

rates(4,8) = R59 * g;
rates(8,4) = R59inv;

rates(4,9) = R56 * n;
rates(9,4) = R56inv;

% node 5
rates(5,9) = R55 * c;
rates(9,5) = R55inv;

rates(5,13) = R54 * n;
rates(13,5) = R54inv;

rates(5,10) = R57 * n;
rates(10,5) = R57inv;

% node 6
rates(6,11) = R37 * c;
rates(11,6) = R37inv;

rates(6,14) = R45 *n;
rates(14,6) = R45inv;

% node 7
rates(7,11) = R38 * g;
rates(11,7) = R38inv;

rates(7,12) = R40 * n;
rates(12,7) = R40inv;

% node 8
rates(8,11) = R39 * n;
rates(11,8) = R39inv;

rates(8,15)  = R46 * n;
rates(15,8) = R46inv;

% node 9
rates(9,12) = R41 * n;
rates(12,9) = R41inv;

rates(9,15) = R47 * g;
rates(15,9) = R47inv;

%node 10
rates(10,14) = R44 * n;
rates(14, 10) = R44inv;

rates(10,15) = R48 * c;
rates(15,10) = R48inv;

% node 11;
rates(11,16) = R33 * n;
rates(16,11) = R33inv;

% node 12
rates(12,13) = R42inv;
rates(13,12) = R42 * c;

rates(12,16) = R35 * g;
rates(16,12) = R35inv;

% node 13
rates(13,14) = R43inv;
rates(14,13) = R43 * g;

% node 14
rates(14,16) = R36 *c;
rates(16,14) = R36inv;

% node 15
rates(15, 16) = R34 * n;
rates(16,15) = R34inv;

rates(2,16) = 210;
rates(16,2) = 70;

%% get k1 and k1inv
%rate_node1_node2 = get_constants_ncoeffs(rate, [node1. node2])
rates_2_8 = get_constants_ncoeffs(rates, [2,8]);
k1 = rates_2_8(2,8);
k1inv = rates_2_8(8,2);

rates_8_15 = get_constants_ncoeffs(rates, [8,15]);
k2 = rates_8_15(8,15);
k2inv = rates_8_15(15,8);

rates_15_16 = get_constants_ncoeffs(rates, [15,16]);
k3 = rates_15_16(15,16);
k3inv = rates_15_16(16,15);

save('ks_GlyT1_top.mat', 'k1', 'k1inv', 'k2', 'k2inv', 'k3', 'k3inv');

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
