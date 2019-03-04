%% MECH 511 - Multigrid Programming Assignment
%% Data Visualisation
% Nicholas Earle

clear;clc;close all;

%% Part 1
% Test function: sin(pi*x)*sin(pi*y)

load('original.dat');       % Original grid
load('OG_coarse.dat');      % Original coarse grid
load('coarse.dat');         % Transferred Coarse grid
load('fine_inj.dat');       % Transferred Fine grid by injection
load('fine_int.dat');       % Transferred Fine grid by interpolation
load('Ecoarse.dat');        % Coarse grid error
load('Einject.dat');        % Injection grid error
load('Einterp.dat');        % Interpolation grid error

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

subplot(2,2,1);
imagesc(xy64, xy64, original);
title('Original Grid');
xlabel('x'); ylabel('y');
axis equal tight

subplot(2,2,2);
imagesc(xy32, xy32, coarse);
title('Coarse Grid');
xlabel('x'); ylabel('y');
axis equal tight

subplot(2,2,3);
imagesc(xy64, xy64, inject);
title('Fine Grid by Injection');
xlabel('x'); ylabel('y');
axis equal tight

subplot(2,2,4);
imagesc(xy64, xy64, interp);
title('Fine Grid by Interpolation');
xlabel('x'); ylabel('y');
axis equal tight

