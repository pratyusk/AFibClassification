clear all; close all; clc;


% load('ecgpuwave_train_data_debug.mat')
% load('Phi_with_sim_P.mat')
load('Phi_train.mat');

features= [2,3,5,7,8,10]; % selected features
labels = Phi(:,1);
Phi = Phi(:,2:end);







%% classify based on 10-fold CV 








test_size = ceil(length(Phi)/10); % 10 for testing each fold
for i = 1:10
    if i*test_size < length(Phi) % might not have exactly even test set size
        vec_test = (i-1)*test_size+1:i*test_size; %1/10 will be test data 
    else
        vec_test = (i-1)*test_size+1:length(Phi);
    end
   % first 10% then next and so on
    vec_train = 1:length(Phi);
    vec_train(vec_test) = []; %rest train data
    %all examples without first test_size than next and so on
    % use linear SVM classifier (function adjusted to get desired outputs!)
    
    
    
    Phi_train_non_PCA_not_normal = Phi(vec_train,:);
    Phi_test_non_PCA_not_normal = Phi(vec_test,:);

    % center/normalize/scale
    Phi_train_non_PCA_not_normal_mean = nanmean(Phi_train_non_PCA_not_normal);  % calculate mean for each column
    Phi_train_non_PCA_not_normal = bsxfun(@minus, Phi_train_non_PCA_not_normal, Phi_train_non_PCA_not_normal_mean);% subtract each column by its mean
    Phi_train_non_PCA_not_normal_std = nanstd(Phi_train_non_PCA_not_normal); % Phi_std is a vector of standard deviation corresponding to each column
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
    
    Phi_PCA_normal = [Phi_train_PCA_normal; Phi_test_PCA_normal];
    labels_PCA = [labels(vec_train); labels(vec_test)];
    %[test_acc_svm_lin_pca_01(i), test_sens_svm_lin_pca_01(i), test_spec_svm_lin_pca_01(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'polynomial',3 , .1,vec_train, vec_test,2);

    %SVM with linear kernel PCA
    [test_acc_svm_lin_pca_01(i), test_sens_svm_lin_pca_01(i), test_spec_svm_lin_pca_01(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'linear', 'auto', .1,vec_train, vec_test,2);
    [test_acc_svm_lin_pca_1(i), test_sens_svm_lin_pca_1(i), test_spec_svm_lin_pca_1(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'linear', 'auto', 1,vec_train, vec_test,2);
    [test_acc_svm_lin_pca_10(i), test_sens_svm_lin_pca_10(i), test_spec_svm_lin_pca_10(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'linear', 'auto', 10,vec_train, vec_test,2);
    [test_acc_svm_lin_pca_100(i), test_sens_svm_lin_pca_100(i), test_spec_svm_lin_pca_100(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'linear', 'auto', 100,vec_train, vec_test,2);
    [test_acc_svm_lin_pca_1000(i), test_sens_svm_lin_pca_1000(i), test_spec_svm_lin_pca_1000(i), ~] = svm_classifier(Phi_PCA_normal(:,features), labels_PCA, 'linear', 'auto', 1000,vec_train, vec_test,2);
    i
    %SVM with linear kernel
    [test_acc_svm_lin_01(i), test_sens_svm_lin_01(i), test_spec_svm_lin_01(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', .1,vec_train, vec_test);
    [test_acc_svm_lin_1(i), test_sens_svm_lin_1(i), test_spec_svm_lin_1(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 1,vec_train, vec_test);
    [test_acc_svm_lin_10(i), test_sens_svm_lin_10(i), test_spec_svm_lin_10(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 10,vec_train, vec_test);
    [test_acc_svm_lin_100(i), test_sens_svm_lin_100(i), test_spec_svm_lin_100(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 100,vec_train, vec_test);
    [test_acc_svm_lin_1000(i), test_sens_svm_lin_1000(i), test_spec_svm_lin_1000(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 1000,vec_train, vec_test);
    i
     %SVM with linear kernel standertized features
    [test_acc_svm_lin_std_01(i), test_sens_svm_lin_std_01(i), test_spec_svm_lin_std_01(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', .1,vec_train, vec_test,1);
    [test_acc_svm_lin_std_1(i), test_sens_svm_lin_std_1(i), test_spec_svm_lin_std_1(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 1,vec_train, vec_test,1);
    [test_acc_svm_lin_std_10(i), test_sens_svm_lin_std_10(i), test_spec_svm_lin_std_10(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 10,vec_train, vec_test,1);
    [test_acc_svm_lin_std_100(i), test_sens_svm_lin_std_100(i), test_spec_svm_lin_std_100(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 100,vec_train, vec_test,1);
    [test_acc_svm_lin_std_1000(i), test_sens_svm_lin_std_1000(i), test_spec_svm_lin_std_1000(i), ~] = svm_classifier(Phi(:,features), labels, 'linear', 'auto', 1000,vec_train, vec_test,1);
    
     
    i
    
    % decision tree
    [test_acc_tree(i), test_sens_tree(i), test_spec_tree(i), ~] = CART_classifier(Phi(:,features), labels,vec_train, vec_test);
    
    %LDA
    [test_acc_lda(i), test_sens_lda(i), test_spec_lda(i), ~] = CART_classifier(Phi(:,features), labels,vec_train, vec_test);
    
    i
end

%% claculate avg and
%std of perfomance measures for each classifier first row mean, second std first 3 svm then RF then CART
Eval_classifier= [ mean(test_acc_svm_lin_01) mean(test_sens_svm_lin_01) mean(test_spec_svm_lin_01);...
    mean(test_acc_svm_lin_1) mean(test_sens_svm_lin_1) mean(test_spec_svm_lin_1);...
    mean(test_acc_svm_lin_10) mean(test_sens_svm_lin_10) mean(test_spec_svm_lin_10);...
    mean(test_acc_svm_lin_100) mean(test_sens_svm_lin_100) mean(test_spec_svm_lin_100);...
    mean(test_acc_svm_lin_1000) mean(test_sens_svm_lin_1000) mean(test_spec_svm_lin_1000);...
    mean(test_acc_svm_lin_std_01) mean(test_sens_svm_lin_std_01) mean(test_spec_svm_lin_std_01);...
    mean(test_acc_svm_lin_std_1) mean(test_sens_svm_lin_std_1) mean(test_spec_svm_lin_std_1);...
    mean(test_acc_svm_lin_std_10) mean(test_sens_svm_lin_std_10) mean(test_spec_svm_lin_std_10);...
    mean(test_acc_svm_lin_std_100) mean(test_sens_svm_lin_std_100) mean(test_spec_svm_lin_std_100);...
    mean(test_acc_svm_lin_std_1000) mean(test_sens_svm_lin_std_1000) mean(test_spec_svm_lin_std_1000);...
    mean(test_acc_svm_lin_pca_01) mean(test_sens_svm_lin_pca_01) mean(test_spec_svm_lin_pca_01);...
    mean(test_acc_svm_lin_pca_1) mean(test_sens_svm_lin_pca_1) mean(test_spec_svm_lin_pca_1);...
    mean(test_acc_svm_lin_pca_10) mean(test_sens_svm_lin_pca_10) mean(test_spec_svm_lin_pca_10);...
    mean(test_acc_svm_lin_pca_100) mean(test_sens_svm_lin_pca_100) mean(test_spec_svm_lin_pca_100);...
    mean(test_acc_svm_lin_pca_1000) mean(test_sens_svm_lin_pca_1000) mean(test_spec_svm_lin_pca_1000);...
     mean(test_acc_tree) mean(test_sens_tree) mean(test_spec_tree);...
      mean(test_acc_lda) mean(test_sens_lda) mean(test_spec_lda)];