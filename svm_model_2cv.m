%% nested cross validation to train svm classifier 

%clear all; close all; clc;
    no_folds_val = 10;

    %rng(1000);
    % for SVM
    Cs = [0.01 0.05 0.1 0.5 1 5 10 20 50 100 1000];
    Sigmas = [0];% 0.01 0.02 0.05 0.1 0.2 0.3 0.5 0.6 0.8 1];

    
    load('Phi_train_new.mat');  
    features = [2 3 5 7 8 10];
    data = Phi(:,[features+1, 1]);
    data(find(data(:,end) == 65),end) = 1;
    data(find(data(:,end) == 78),end) = 0;
    
%     ratio = sum(data(:,end))/(size(data,1)-sum(data(:,end)));
%     if(ratio>=1)        
%         data = [data ; repmat(data(data(:,end)==0,:),round(ratio-1),1)];
%     else
%         data = [data ; repmat(data(data(:,end)==1,:),round(1/ratio-1),1)];
%     end        
    
    data0 = data(data(:,end)==0,:);
    data1 = data(data(:,end)==1,:);
    
    rperm0 = randperm(length(data0));
    rperm1 = randperm(length(data1));
    data0 = data0(rperm0,:);
    data1 = data1(rperm1,:);
    
    fold_size0 = ceil(size(data0,1)/no_folds_val);
    fold_size1 = ceil(size(data1,1)/no_folds_val);       

    accs = zeros(1,no_folds_val);
    TPs = zeros(1,no_folds_val);
    TNs = zeros(1,no_folds_val);
    FPs = zeros(1,no_folds_val);
    FNs = zeros(1,no_folds_val);
    lengths = zeros(1,no_folds_val);
    best_cs = zeros(1,no_folds_val);
    best_ss = zeros(1,no_folds_val);
    p_labels = cell(1,no_folds_val);%?
    p_scores = cell(1,no_folds_val);%?
    a_labels = cell(1,no_folds_val);%?
    parfor fi = 1:no_folds_val
        tc = tic;
        
        fold_test0 = data0(fold_size0*(fi-1)+1:min(fold_size0*fi,size(data0,1)),:);
        fold_test1 = data1(fold_size1*(fi-1)+1:min(fold_size1*fi,size(data1,1)),:);

        fold_train0 = data0([1:fold_size0*(fi-1),min(fold_size0*fi,size(data0,1))+1:end],:);
        fold_train1 = data1([1:fold_size1*(fi-1),min(fold_size1*fi,size(data1,1))+1:end],:);            

        fold_test = [fold_test0 ; fold_test1];
        fold_test = fold_test(randperm(size(fold_test,1)),:);

        fold_train = [fold_train0 ; fold_train1];
        fold_train = fold_train(randperm(size(fold_train,1)),:);
        
        x_train = fold_train(:,1:end-1);
        y_train = fold_train(:,end);
        x_test = fold_test(:,1:end-1);
        y_test = fold_test(:,end);
        
        acc_fold = zeros(length(Cs),length(Sigmas));
        w = 1;%*sum(y_train==0)/sum(y_train==1);
        for ci = 1:length(Cs)
            for si = 1:length(Sigmas)                         

                if(Sigmas(si) == 0) 
                    svm_model = fitcsvm(x_train, y_train, 'KernelFunction', 'linear', 'PolynomialOrder', [], 'KernelScale', 'auto', 'BoxConstraint', Cs(ci), 'Standardize', 1, 'Crossval','on', 'Cost', [0, 1 ; w, 0]);
                else
                    svm_model = fitcsvm(x_train, y_train, 'KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', Sigmas(si), 'BoxConstraint', Cs(ci), 'Standardize', 1, 'Crossval','on', 'Cost', [0, 1 ; w, 0]);
                end                                      
                
                acc_fold(ci,si) = svm_model.kfoldLoss();              
            end
        end
        [mn, mi] = min(acc_fold);
        [~, best_si] = min(mn);
        best_ci = mi(best_si);
        
        if(Sigmas(best_si) == 0)
            svm_model = fitcsvm(x_train, y_train, 'KernelFunction', 'linear', 'PolynomialOrder', [], 'KernelScale', 'auto', 'BoxConstraint', Cs(best_ci), 'Standardize', 1, 'Cost', [0, 1 ; w, 0]);
        else
            svm_model = fitcsvm(x_train, y_train, 'KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', Sigmas(best_si), 'BoxConstraint', Cs(best_ci), 'Standardize', 1, 'Cost', [0, 1 ; w, 0]);
        end 
        
        [labels, scores] = predict(svm_model, x_test);
        p_labels{fi} = labels;
        p_scores{fi} = scores;
        a_labels{fi} = y_test;
        best_cs(fi) = best_ci;
        best_ss(fi) = best_si;
        
        
        lengths(fi) = length(p_labels{fi});
        accs(fi) = 100*(1-sum(abs(p_labels{fi}-a_labels{fi}))./length(p_labels{fi}));
        TPs(fi) = sum(a_labels{fi}>0 & p_labels{fi}>0);
        TNs(fi) = sum(a_labels{fi}==0 & p_labels{fi}==0);
        FPs(fi) = sum(a_labels{fi}==0 & p_labels{fi}>0);
        FNs(fi) = sum(a_labels{fi}>0 & p_labels{fi}==0);    
        toc(tc);
    end
    
    pred_labels = [];
    for fi = 1:no_folds_val
        pred_labels = [pred_labels; p_labels{fi}];
    end
    
    disp(accs)
    disp(sum(accs.*lengths)/sum(lengths))
    disp(std(accs))
    
    save('resutls','p_labels','p_scores','a_labels','best_cs','best_ss');
