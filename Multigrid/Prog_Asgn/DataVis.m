%% MECH 511 - Multigrid Programming Assignment
%% Data Visualisation
% Nicholas Earle

clear;clc;close all;

%% Part 1
% Test function: sin(pi*x)*sin(pi*y)

load('Data/original.dat');       % Original grid
load('Data/OG_coarse.dat');      % Original coarse grid
load('Data/coarse.dat');         % Transferred Coarse grid
load('Data/fine_inj.dat');       % Transferred Fine grid by injection
load('Data/fine_int.dat');       % Transferred Fine grid by interpolation
load('Data/Ecoarse.dat');        % Coarse grid error
load('Data/Einject.dat');        % Injection grid error
load('Data/Einterp.dat');        % Interpolation grid error

original = original(2:end-1, 2:end-1);
OG_coarse = OG_coarse(2:end-1, 2:end-1);
coarse = coarse(2:end-1, 2:end-1);
inject = fine_inj(2:end-1, 2:end-1);
interp = fine_int(2:end-1, 2:end-1);
Ecoarse = Ecoarse(2:end-1, 2:end-1);
Einject = Einject(2:end-1, 2:end-1);
Einterp = Einterp(2:end-1, 2:end-1);

xy32 = linspace(0,1,32);
xy64 = linspace(0,1,64);

figure(1);
imagesc(xy64, xy64, original);
title('Original Grid');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

figure(2);
imagesc(xy32, xy32, coarse);
title('Coarse Grid');
xlabel('x'); ylabel('y');
axis equal tight
colorbar

figure(3);
imagesc(xy64, xy64, inject);
title('Fine Grid by Injection');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

figure(4);
imagesc(xy64, xy64, interp);
title('Fine Grid by Interpolation');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

figure(5);
imagesc(xy32, xy32, Ecoarse);
title('Error in Coarse Grid Transfer, L2 = 3.0E-4');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

figure(6);
imagesc(xy64, xy64, Einject);
title('Error in Injection Grid Transfer, L2 = 1.7E-2');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

figure(7);
imagesc(xy64, xy64, Einterp);
title('Error in Interpolation Grid Transfer, L2 = 1.2E-3');
xlabel('x'); ylabel('y');
axis equal tight
colorbar;

%% Part 3

load('Data/L2_1.dat'); L2_1 = L2_1(2:end);
load('Data/L2_2.dat'); L2_2 = L2_2(2:end);
load('Data/L2_3.dat'); L2_3 = L2_3(2:end);
load('Data/L2_4.dat'); L2_4 = L2_4(2:end);
load('Data/L2_5.dat'); L2_5 = L2_5(2:end);

figure();
semilogy(L2_1, 'LineWidth', 1.5);
hold on;
semilogy(L2_2, 'LineWidth', 1.5);
semilogy(L2_3, 'LineWidth', 1.5);
semilogy(L2_4, 'LineWidth', 1.5);
semilogy(L2_5, 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('L_2 Norm');
title('Convergence Rate for simple V-Cycles');
ylim([1E-9 1E0]);
legend('1 Mesh - 8314','2 Meshes - 1386','3 Meshes - 321','4 Meshes - 79','5 Meshes - 20');

%% Part 4

load('Data/L2_w0.500000.dat'); L2_w05 = L2_w0_500000(2:end);
load('Data/L2_w0.600000.dat'); L2_w06 = L2_w0_600000(2:end);
load('Data/L2_w0.666667.dat'); L2_w067 = L2_w0_666667(2:end);
load('Data/L2_w0.800000.dat'); L2_w08 = L2_w0_800000(2:end);
load('Data/L2_w1.000000.dat'); L2_w10 = L2_w1_000000(2:end);
load('Data/L2_w1.250000.dat'); L2_w125 = L2_w1_250000(2:end);
load('Data/L2_w1.500000.dat'); L2_w15 = L2_w1_500000(2:end);
load('Data/L2_double.dat'); L2_double = L2_double(2:end);

figure();
semilogy(L2_w05, 'LineWidth', 1.5);
hold on;
semilogy(L2_w06, 'LineWidth', 1.5);
semilogy(L2_w067, 'LineWidth', 1.5);
semilogy(L2_w08, 'LineWidth', 1.5);
semilogy(L2_w10, 'LineWidth', 1.5);
semilogy(L2_w125, 'LineWidth', 1.5);
semilogy(L2_w15, 'LineWidth', 1.5);
semilogy(L2_double, 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('L_2 Norm');
title('Convergence Rate for a 4-level V-cycle with differing \omega');
ylim([1E-9 1E0]);
legend('\omega = 0.5 - 108','\omega = 0.6 - 89','\omega = 0.67 - 79',...
    '\omega = 0.8 - 62','\omega = 1.0 - 44','\omega = 1.25 - 26',...
    '\omega = 1.5 - 10', 'Double Pass - 36x2=72');



