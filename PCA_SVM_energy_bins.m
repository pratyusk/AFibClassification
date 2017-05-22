%% PCA for feature selection

clear all;
close all;
clc;
%load data
load('Phi_train.mat');
train_size = round(size(Phi, 1)/2);
% train_size = 900;
% extract features 
labels = Phi(:,1); % extract all labels
labels_train = labels(1:train_size); % training labels
labels_test = labels(train_size+1:end); % testing labels

Phi = Phi(:,2:end); % features without the labels
% Phi = Phi(:,[11:15,18]);
Phi_non_PCA_not_normal = Phi;
Phi_train_non_PCA_not_normal = Phi_non_PCA_not_normal(1:train_size,:);
Phi_test_non_PCA_not_normal = Phi_non_PCA_not_normal(train_size+1:end,:);

% center/normalize/scale
Phi_train_non_PCA_not_normal_mean = mean(Phi_train_non_PCA_not_normal);  % calculate mean for each column
Phi_train_non_PCA_not_normal = bsxfun(@minus, Phi_train_non_PCA_not_normal, Phi_train_non_PCA_not_normal_mean);% subtract each column by its mean
Phi_train_non_PCA_not_normal_std = std(Phi_train_non_PCA_not_normal); % Phi_std is a vector of standard deviation corresponding to each column
Phi_train_non_PCA_not_normal = bsxfun(@rdivide, Phi_train_non_PCA_not_normal, Phi_train_non_PCA_not_normal_std); % divide each column element in Phi by Phi_std (corresponding index)
Phi_train_non_PCA_not_normal(isinf(Phi_train_non_PCA_not_normal)) = 0; % set all divide by 0 to 0
Phi_train_non_PCA_not_normal(isnan(Phi_train_non_PCA_not_normal)) = 0; %

% center/normalize/scale
Phi_test_non_PCA_not_normal = bsxfun(@minus, Phi_test_non_PCA_not_normal, Phi_train_non_PCA_not_normal_mean); % subtract each column by its mean
Phi_test_non_PCA_not_normal = bsxfun(@rdivide, Phi_test_non_PCA_not_normal, Phi_train_non_PCA_not_normal_std); % divide each column element in Phi by Phi_std (corresponding index)
Phi_test_non_PCA_not_normal(isinf(Phi_test_non_PCA_not_normal)) = 0; % set all divide by 0 to 0
Phi_test_non_PCA_not_normal(isnan(Phi_test_non_PCA_not_normal)) = 0; %

% split data Phi_train, Phi_test
Phi_train_non_PCA_normal = Phi_train_non_PCA_not_normal;
Phi_test_non_PCA_normal = Phi_test_non_PCA_not_normal;
% Phi_train = Phi(1:304,:);
% Phi_test = Phi(304:608,:);

% apply PCA
[coeff,score,latent] = pca(Phi_train_non_PCA_normal);
Phi_train_PCA_normal = Phi_train_non_PCA_normal * coeff;
Phi_test_PCA_normal = Phi_test_non_PCA_normal * coeff;

% Phi_train_PCA_normal = Phi_train_non_PCA_normal * 2;
% Phi_test_PCA_normal = Phi_test_non_PCA_normal * 2;

% test with SVM
[test_acc_PCA, train_acc_PCA, SVM_Model_PCA_normal] = svm_classifier(Phi_train_PCA_normal,labels_train,'polynomial',1,1, Phi_test_PCA_normal, labels_test);
[test_acc_normal, train_acc_normal, SVM_Model_normal] = svm_classifier(Phi_train_non_PCA_normal,labels_train,'polynomial',1,1, Phi_test_non_PCA_normal, labels_test);
[test_acc_not_normal, train_acc_not_normal, SVM_Model_not_normal] = svm_classifier(Phi_train_non_PCA_not_normal,labels_train,'polynomial',1,1,Phi_test_non_PCA_not_normal, labels_test);
% figure();
% svm_plot_classes(Phi_train_PCA_normal,labels_train,SVM_Model_PCA_normal);
% Phi_pca_train + Phi_train *  coeff

%svm(Phi_pca) --> predict coeff*Phi_test, 
%svm(Phi) --> predict Phi_test
