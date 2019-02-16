%% MECH 511 - HW2 - Unstructured Mesh
% Nick Earle

clear; close all; clc;

%%

P = [3.25 3 2 4.5 1.5 2.5; 1.25 0.25 1 0.5 0 2.25];
T = [1 1 2 1; 2 2 3 3; 3 4 5 6];

Temp = [100 102 97 101];

C = [mean(P(1,T(:,1))) mean(P(1,T(:,2))) mean(P(1,T(:,3))) mean(P(1,T(:,4)))
    mean(P(2,T(:,1))) mean(P(2,T(:,2))) mean(P(2,T(:,3))) mean(P(2,T(:,4)))];

D = zeros(3, 3);
for i = 2:4
    D(1, i-1) = C(1,i) - C(1,1);
    D(2, i-1) = C(2,i) - C(2,1);
    D(3, i-1) = Temp(i) - Temp(1);
end

A = D(1:2, :)';
b = D(end, :)';

dT = (A'*A)^-1 * (A'*b);
