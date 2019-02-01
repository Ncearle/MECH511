%% MECH 511 - Multigrid Assignment

clear; clc; close all;

xmin = 0; xmax = 1;
imax = 1000;
dx = (xmax - xmin) / imax;

k = imax/2:imax;
sig = @(w) max(abs(1 - w*(1 - cos(k*pi*dx))));


sig2 = @(w1) (1 - w1*(1 - cos(k*pi*dx))) .* (1 - w1/(3*w1-1)*(1 - cos(k*pi*dx)));
g = @(w) abs(max(sig2(w)) - abs(min(sig2(w))));

x = fminsearch(g, 0.5);

% sig3 = @(w1, w2) (1 - w1*(1 - cos(imax/2*pi*dx))) .* (1 - w2*(1 - cos(imax/2*pi*dx)))...
%     - (1 - w1*(1 - cos(imax*pi*dx))) .* (1 - w2*(1 - cos(imax*pi*dx)));


sig2 = @(w1, w2) (1 - w1*(1 - cos(k*pi*dx))) .* (1 - w2*(1 - cos(k*pi*dx)));

w1 = x;
w2 = w1./(3*w1-1);



plot(k, sig2(w1, w2));
hold on

legend(sprintf('w1=%0.4f, w2=%0.4f', w1, w2));


% y = fminsearch(sig3, w);
% 
% x = fminsearch(sig, 0.6);
% 
% plot(k, sig2(y));
% hold on
% plot(k, sig2(0.6, 0.8, k));
% plot(k, sig2(0.6, 1, k));
% plot(k, sig2(0.55, 0.85, k));

% figure(1);
% for w = 1:10
%     plot(k, sig(w/10));
%     hold on;
%     q = w/10;
%     leg{w} = sprintf('w = %0.1f',q);
% end
% legend(leg);

