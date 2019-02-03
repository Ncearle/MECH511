%% MECH 511 - Multigrid Assignment

clear; clc; close all;

xmin = 0; xmax = 1;
imax = 100;
dx = (xmax - xmin) / imax;

k = imax/2:imax;
sig = @(w) 1 - w*(1 - cos(k*pi*dx));


sig2 = @(w1) (1 - w1*(1 - cos(k*pi*dx))) .* (1 - w1/(3*w1-1)*(1 - cos(k*pi*dx)));
g = @(w) abs(max(sig2(w)) - abs(min(sig2(w))));

x = fminsearch(g, 0.5);

% sig3 = @(w1, w2) (1 - w1*(1 - cos(imax/2*pi*dx))) .* (1 - w2*(1 - cos(imax/2*pi*dx)))...
%     - (1 - w1*(1 - cos(imax*pi*dx))) .* (1 - w2*(1 - cos(imax*pi*dx)));


sig2 = @(w1, w2) (1 - w1*(1 - cos(k*pi*dx))) .* (1 - w2*(1 - cos(k*pi*dx)));
% 
w1 = x;
w2 = w1./(3*w1-1);
% 
% w1 = 0.6; w2 = 0.8;
sig2 = sig2(w1, w2);
sig = sig(2/3).^2;
plot(k, sig2, k, sig);
title('Amplification factor vs. wave number');
xlabel('Wave number');
ylabel('Amplification factor');
% 
legend(sprintf('\\omega_{1}=%0.4f, \\omega_{2}=%0.4f', w1, w2), '\omega = 0.6667');

% x = fminsearch(sig, 0.6);

% figure(1);
% for w = 5:10
%     plot(k, sig(w/10));
%     hold on;
%     q = w/10;
%     leg{w-4} = sprintf('\\omega = %0.1f',q);
% end
% legend(leg);
% title('Amplification factor vs. wave number');
% xlabel('Wave number');
% ylabel('Amplification factor');
% 
% sig = @(w) 1 - w*(1 - cos(k*pi*dx));
% plot(k, sig(x));
% legend(sprintf('\\omega = %0.4f', x));
% title('Optimal Amplification factor vs. wave number');
% xlabel('Wave number');
% ylabel('Amplification factor');

