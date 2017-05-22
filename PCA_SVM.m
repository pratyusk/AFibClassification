%% PCA for feature selection

clear all;
close all;
clc;
%load data
load('Phi_2.mat');
train_size = round(size(Phi, 1)/2);
% extract features 
labels = Phi(:,1); % extract all labels
labels_train = labels(1:train_size); % training labels
labels_test = labels(train_size+1:end); % testing labels

Phi = Phi(:,2:end); % features without the labels
Phi = Phi(:,[1:5,11:19]);
Phi_non_PCA_not_normal = Phi;
Phi_train_non_PCA_not_normal = Phi_non_PCA_not_normal(1:train_size,:);
Phi_test_non_PCA_not_normal = Phi_non_PCA_not_normal(train_size+1:end,:);

% center/normalize/scale
Phi = detrend(Phi, 'constant'); % subtract each column by its mean
Phi_std = std(Phi); % Phi_std is a vector of standard deviation corresponding to each column
Phi = bsxfun(@rdivide, Phi, Phi_std); % divide each column element in Phi by Phi_std (corresponding index)
Phi(isinf(Phi)) = 0; % set all divide by 0 to 0
Phi(isnan(Phi)) = 0; %

% split data Phi_train, Phi_test
Phi_train_non_PCA_normal = Phi(1:train_size,:);
Phi_test_non_PCA_normal = Phi(train_size+1:end,:);
% Phi_train = Phi(1:304,:);
% Phi_test = Phi(304:608,:);

% apply PCA
[coeff,score,latent] = pca(Phi_train_non_PCA_normal);
Phi_train_PCA_normal = Phi_train_non_PCA_normal * coeff';
Phi_test_PCA_normal = Phi_test_non_PCA_normal * coeff';

% test with SVM
[test_acc_PCA, train_acc_PCA, SVM_Model_PCA_normal] = svm_classifier(Phi_train_PCA_normal,labels_train,'linear',1,100, Phi_test_PCA_normal, labels_test);
[test_acc_normal, train_acc_normal, SVM_Model_normal] = svm_classifier(Phi_train_non_PCA_normal,labels_train,'linear',1,100, Phi_test_non_PCA_normal, labels_test);
[test_acc_not_normal, train_acc_not_normal, SVM_Model_not_normal] = svm_classifier(Phi_train_non_PCA_not_normal,labels_train,'linear',1,100,Phi_test_non_PCA_not_normal, labels_test);
% figure();
% svm_plot_classes(Phi_train_PCA_normal,labels_train,SVM_Model_PCA_normal);
% Phi_pca_train + Phi_train *  coeff

%svm(Phi_pca) --> predict coeff*Phi_test, 
%svm(Phi) --> predict Phi_test
