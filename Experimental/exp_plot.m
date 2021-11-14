
[sheet1, headings1] = xlsread("GlyR-GlyT data.xlsx" ,1);
[sheet2, headings2] = xlsread("GlyR-GlyT data.xlsx" ,2);


% sodium dependence of glycine at 10\muM
times = sheet2(:,end-3);
Na96= sheet2(:,end-2);
Na10= sheet2(:,end-1);
Na1= sheet2(:,end);

f = figure(1); hold on;
% plot(times, Na96, times, Na10, times, Na1)
% xlabel("time(s)")
% ylabel("Current (nA)")
% title("Current in stop GlyR/GlyT2 flow experiment at different sodium concentrations")
% l = legend(["96mM", "10mM", "1mM"]);
% title(l,"Extracellular Sodium");


%% isolate sections


stop96_i = 403;
stop96_f = 488;

t96s = times(stop96_i:stop96_f);


Na96s = Na96(stop96_i:stop96_f);
Na96d = Na96s -min(Na96s);
Na96d = Na96d / max(Na96d);
%Na96d = Na96d / max(Na96d);

H96 = Na96d ./(1.02- Na96d);



stop10_i = 901;
stop10_f = 1009;

t10s = times(stop10_i:stop10_f);


Na10s = Na10(stop10_i:stop10_f);
Na10d = Na10s -min(Na10s);
Na10d = Na10d / max(Na10d);

H10 = Na10d ./(1.02- Na10d);



stop1_i = 1378;
stop1_f = 1479;

t1s = times(stop1_i:stop1_f);

Na1s = Na1(stop1_i:stop1_f);
Na1d = Na1s -min(Na1s);
Na1d = Na1d / max(Na1d);

H1 = Na1d ./(1.02- Na1d);



% start time fromjust after zero (dont want to break log(x))
t96s = t96s(2:end) - t96s(1); 
t1s = t1s(2:end) - t1s(1);
t10s = t10s(2:end) - t10s(1);

H96 = H96(2:end);
H10 = H10(2:end);
H1 = H1(2:end);

f1 = fit(log(t96s), log(H96), 'poly1');
f2 = fit(log(t10s), log(H10), 'poly1');
f3 = fit(log(t1s), log(H1), 'poly1');

figure(3)
plot(t96s / max(t96s), Na96s(2:end) / min(Na96s), t10s / max(t10s), Na10s(2:end) / min(Na10s), t1s  / max(t1s), Na1s(2:end) / min(Na1s))
title("plot of normalised membrane current decay")

l = legend(["96", "10", "1"]);
title(l,"Extracellular Sodium (mM)");
xlabel("t / t_{max}")
ylabel("I / I_{max}")


% figure(2)
% loglog(t96s, H96, t10s, H10, t1s, H1)
% title("log scale plot of normalised membrane current decay")
% 
% l = legend(["96mM, q = (2.089, 2.15)", "10mM, q = (1.918, 1.934)", "1mM, q = (2.124, 2.157)"]);
% title(l,"Extracellular Sodium");
% xlabel("time (s)")
% ylabel("J / (J_{max} - J)")

%% glycine measurement

% load('Glye150.mat') % this is from the homestatic conditions run
% figure(3)
% cutoff = 40;
% t = t(cutoff:end);
% GLY_e = GLY_e(cutoff:end);
% GLY_e = GLY_e / max(GLY_e);
% 
% HGe = GLY_e ./ (1.06 - GLY_e);
% loglog(t, HGe);
% fge = fit(log(t), log(HGe), 'poly1');
% xlabel("dimensionless time")
% ylabel("G_e / (G_{max} -G), G dimensionless")
% title("Log-Log plot of transporter under typical conditions")
