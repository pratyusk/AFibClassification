%% train model with best chosen parameter from nested cv 'svm_model_2cv.m'
%train model on all data for submission
clear all;
close all;
clc;
load('Phi_train_new.mat');
Phi_train = Phi;
load('Phi_test_new.mat');
Phi_train = [Phi_train;Phi];


C = 0.05;
w = 1;
features = [2 3 5 7 8 10];
x_train = Phi_train(:,features+1);
y_train = Phi_train(:,1) == 65;



%% train model on 100% train set
svm_model = fitcsvm(x_train, y_train, 'KernelFunction', 'linear', 'PolynomialOrder', [], 'KernelScale', 'auto', 'BoxConstraint', C, 'Standardize', 1, 'Cost', [0, 1 ; w, 0]);






