%% train model with best chosen parameter from nested cv 'svm_model_2cv.m'
%train model on 70% train set and 30% test set
clear all;
close all;
clc;
load('Phi_train_new.mat');
Phi_train = Phi;
load('Phi_test_new.mat');
Phi_test = Phi;

delta = 0.3428;
gamma = 0.4900;
features = [2 3 5 7 8 10];
x_train = Phi_train(:,features+1);
y_train = Phi_train(:,1) == 65;

x_test = Phi_test(:,features+1);
y_test = Phi_test(:,1) == 65;


%% Applying PCA
% Phi = detrend(Phi, 'constant'); % subtract each column by its mean
x_mean = nanmean(x_train);
x_train = bsxfun(@minus, x_train, x_mean);
x_train_std = nanstd(x_train); % Phi_std is a vector of standard deviation corresponding to each column
x_train = bsxfun(@rdivide, x_train, x_train_std); % divide each column element in Phi by Phi_std (corresponding index)
x_train(isinf(x_train)) = 0; % set all divide by 0 to 0
x_train(isnan(x_train)) = 0; %

x_test = bsxfun(@minus, x_test, x_mean);
x_test = bsxfun(@rdivide, x_test, x_train_std); % divide each column element in Phi by Phi_std (corresponding index)
x_train(isinf(x_train)) = 0; % set all divide by 0 to 0
x_train(isnan(x_train)) = 0; %

% apply PCA
[coeff,score,latent] = pca(x_train);
x_train_pca = x_train * coeff';
x_test_pca = x_test * coeff';




%% train model on 70% train set
lda_model = fitcdiscr(x_train_pca,y_train,'Delta',delta,'Gamma',gamma);



%% test model
for i = 1:length(y_test)
   if length(find(isnan(x_test_pca(i,:)))) > 0
       test_pred_y(i) = 1;
       continue
   end
   test_pred_y(i) = predict(lda_model, x_test_pca(i,:));       
end
%% result
y_test = y_test';
%true pos is class 1 classified as class 1
 true_pos_test = sum(y_test == test_pred_y & y_test == 1);


% true negative is class 0 classified as class 0
true_neg_test = sum(y_test == test_pred_y & y_test == 0);


% false positive is class 0 classified as class 1
false_pos_test = sum(y_test ~= test_pred_y &  y_test == 1);


% false negative is class 1 classified as class 0
false_neg_test = sum(y_test ~= test_pred_y &  y_test == 0);
    
% accuracy = (true pos + true neg)/(true pos + true neg + false pos + false neg)
test_acc = (true_pos_test + true_neg_test)/(true_pos_test + true_neg_test + false_pos_test + false_neg_test)

% sensitivity = true pos/(true pos + false neg) = true pos/ pos class
test_sens = true_pos_test/(true_pos_test + false_neg_test)

% Specificity = true neg/(true neg + false pos)
test_spec = true_neg_test/(true_neg_test + false_pos_test)




