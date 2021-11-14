function [k1, k2,k3,k4, kinv1, kinv2, kinv3, kinv4] =  GlyT2_top(n,c,g)

% define Rs here, node numbering is done in GLYT2_top.png
num_nodes = 21;
coop = 8;
R49 = 1e3;
R49inv = 300 * coop;

R64inv = 100* coop;
R64 = 1e6;

R60 = 1e6;
R60inv = 1e4 * coop;

 R58 = 1e4;
R58inv = 100 * coop;

% node 2
R61 = 1e3;
R61inv = 300*coop;

 R62 =1e6;
R62inv = 1e4 * coop;

R63 = 1e4;
R63inv = 100* coop;

% node 3
 R50 = 1e6;
 R50inv = 100* coop;

 R52 = 1e6;
 R52inv = 1e4*coop;

R53 = 1e4;
R53inv = 100*coop;

%node 4
R51 =  1e3;
R51inv = 300*coop;

R59 = 10^6;
R59inv = 100 * coop;

R56 = 10^4;
R56inv = 100*coop;

% node 5
R54 = 10^3;
R54inv = 300*coop;

R55 =1e6;
R55inv = 1e4 * coop;

R57 = 10^6;
R57inv = 100* coop;

% node 6
R37 =10^6;
R37inv = 1^4 * coop;

R45 = 10^4;
R45inv = 100*coop;

% node 7
R38 = 10^6;
R38inv = 100*coop;

R40 = 10^4;
R40inv = 100*coop;

% node 8
R39 = 10^3;
R39inv = 300*coop;

R46 = 10^6;
R46inv = 10^4 * coop;

% node 9
R41 = 10^3;
R41inv = 300*coop;

R47 = 10^6;
R47inv = 100*coop;

%node 10
R44 = 10^3;
R44inv = 300*coop;

R48 = 10^6;
R48inv = 10^4 * coop;

% node 11;
R79 = 10^4;
R79inv = 100*coop;

R33 = 10^4;
R33inv = 100*coop;

% node 12
R35 = 10^6;
R35inv = 100*coop;

R81 = 10^4;
R81inv = 100*coop;

R42inv = 10^4 * coop;
R42 = 10^6;

% node 13
R43inv = 100*coop;
R43 = 10^6;

% node 14
R36 =10^6;
R36inv = 10^4 * coop;

R82 = 10^4;
R82inv = 100*coop;

% node 15
R80 = 10^4;
R80inv = 100*coop;

R34 = 10^3;
R34inv = 300*coop;

% node 16
R72 =10^4;
R72inv = 100;

%node 17
R73 = 10^3;
R73inv = 300;

%node 18
R74 = 10^4;
R74inv = 100;

% node 19
R71 =10^6;
R71inv = 10^4;

%node 20
R70 =10^6;
R70inv = 100;
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
rates(5,13) = R54 * n;
rates(13,5) = R54inv;

rates(5,9) = R55 * c;
rates(9,5) = R55inv;

rates(5,10) = R57 * g;
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
rates(11,16) = R79 *n;
rates(16,11) = R79inv;

rates(11,18) = R33 * n;
rates(18,11) = R33inv;

% node 12
rates(12,18) = R35 * g;
rates(18,12) = R35inv;

rates(12,20) = R81 * n;
rates(20,12) = R81inv;

rates(12,13) = R42inv;
rates(13,12) = R42 * c;

% node 13
rates(13,14) = R43 *g;
rates(14,13) = R43inv;

% node 14
rates(14,18) = R36 *c;
rates(18,14) = R36inv;

rates(14,19) = R82 * n;
rates(19,14) = R82inv;

% node 15
rates(15,17) = R80 * n;
rates(17,15) = R80inv;

rates(15, 18) = R34 * n;
rates(18,15) = R34inv;

% node 16
rates(16, 21) = R72 * n;
rates(21, 16) = R72inv;

%node 17
rates(17,21) = R73 * n;
rates(21,17) = R73inv;

%node 18
rates(18,21) = R74 * n;
rates(21,18) = R74inv;

% node 19
rates(19,21) = R71 * c;
rates(21,19) = R71inv;

%node 20
rates(20,21) = R70 * g;
rates(21,20) = R70inv;

%% Extra pathway for pseudo-equilibrium calculation
rates(2,21) = 300;
rates(21,2) = 300;


%% get k and kinv values
rates_2_8 = get_constants_ncoeffs(rates, [2,  8]);

k1 = rates_2_8(2,8);
kinv1 = rates_2_8(8,2);

rates_8_15 = get_constants_ncoeffs(rates, [8, 15]);
k2 = rates_8_15(8,15);
kinv2 = rates_8_15(15,8);

rates_15_16 = get_constants_ncoeffs(rates, [15, 16]);
k3 = rates_15_16(15,16);
kinv3 = rates_15_16(16,15);

rates_16_21 = get_constants_ncoeffs(rates, [16,  num_nodes]);
k4 = rates_16_21(16, 21);
kinv4 = rates_16_21(21,16);



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
