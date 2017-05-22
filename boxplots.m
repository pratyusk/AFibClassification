clear all;
close all; clc;
load('Phi_train_new.mat');
load('feature_labels.mat');
labels = Phi(:,1);
Phi = Phi(:,2:end);
for i=1:size(Phi,2) % loop till number of features
    figure()
    boxplot(Phi(:,i),labels)
    ylabel('Values')
    xlabel('Classes')
    title(strcat('Boxplot of feature_ ',feature_labels(i,:)),'Interpreter', 'none')
%      saveas(gcf,strcat('final Boxplot of feature_ ',feature_labels(i,:)),'png')
%     savefig(strcat('final Boxplot of feature_ ',feature_labels(i,:)));
end