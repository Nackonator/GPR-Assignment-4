%This file is meant for debugging and testing DELETE BEFORE HANDING IN
clc
clear all
close all
addpath("functions")
addpath("test functions")
%% Testing kernel_matrix function
hyp = [0.05;1]; %defining hyper parameters
[X,y,Xstar] = getdata(5,hyp);

K1 = kernel_matrix(X,hyp(1),hyp(2));
K2 = kernel_matrix_unopt(X,hyp(1),hyp(2));

norm(K1-K2);
t = sum(K1-K2 ~= 0);