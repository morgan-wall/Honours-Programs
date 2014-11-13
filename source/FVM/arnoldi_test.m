clear all; 
close all;

mex arnoldi.c

load('arnoldi_test_simple.mat');

restart_value = 5;

% for i = 1:100
    [xresult, i] = ...
        gmres_general(A, b, x0, L, U, max_iter, restart_value, error_tol);
% end