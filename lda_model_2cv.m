%% nested cross validation to train svm classifier 

%clear all; close all; clc;
    no_folds_val = 10;

    %for LDA
    Deltas = [.05 0.1 0.5 1 5 10 20 50 100 1000];
    Gammas = [0 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    
    load('Phi_train_new.mat');  
    features = [2 3 5 7 8 10];
    data = Phi(:,[features+1, 1]);
    data(find(data(:,end) == 65),end) = 1;
    data(find(data(:,end) == 78),end) = 0;
    
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
    best_gs = zeros(1,no_folds_val);
    best_ds = zeros(1,no_folds_val);
    p_labels = cell(1,no_folds_val);%?
    p_scores = cell(1,no_folds_val);%?
    a_labels = cell(1,no_folds_val);%?
    
    accs_qda = zeros(1,no_folds_val);
    TPs_qda = zeros(1,no_folds_val);
    TNs_qda = zeros(1,no_folds_val);
    FPs_qda = zeros(1,no_folds_val);
    FNs_qda = zeros(1,no_folds_val);
    lengths_qda = zeros(1,no_folds_val);
    p_labels_qda = cell(1,no_folds_val);%?
    p_scores_qda = cell(1,no_folds_val);%?
    a_labels_qda = cell(1,no_folds_val);%?
    
    for fi = 1:no_folds_val
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
        
        acc_fold = zeros(length(Deltas),1);
        w = 1;%*sum(y_train==0)/sum(y_train==1);
        
        %% inner CV to determine best Hyperparameters
        lda_model = fitcdiscr(x_train,y_train,'SaveMemory','on','FillCoeffs','off');
        [err,gamma,delta,numpred] = cvshrink(lda_model,'NumGamma',10,'NumDelta',10);
        
%         [mn, mi] = min(acc_fold);
%         [~, best_si] = min(mn);
%         best_ci = mi(best_si);
        minerr = min(min(err));
        [p,q] = find(err == minerr); % Subscripts of err producing minimal error
        idx = sub2ind(size(delta),p,q); % Convert from subscripts to linear indices
%         lowerr = min(min(err));
%         lownum = min(min(numpred));
%         low200 = min(min(err(numpred <= 200)));
%         lownum = min(min(numpred(err == low200)));
%         [r,s] = find((err == lowerr) & (numpred == lownum));
%         [r,s] = find((err == low200) & (numpred == lownum));
%         lda_model.Gamma = gamma(p(1));
%         lda_model.Delta = delta(idx(1)); 

           %% train LDA based on best parameters from inner CV
        lda_model = fitcdiscr(x_train,y_train,'Delta',delta(idx(1)),'Gamma',gamma(p(1)));

        
        
        [labels, scores] = predict(lda_model, x_test);
        p_labels{fi} = labels;
        p_scores{fi} = scores;
        a_labels{fi} = y_test;
        best_gs(fi) = gamma(p(1));
        best_ds(fi) = delta(idx(1));
        
        
        lengths(fi) = length(p_labels{fi});
        accs(fi) = 100*(1-sum(abs(p_labels{fi}-a_labels{fi}))./length(p_labels{fi}));
        TPs(fi) = sum(a_labels{fi}>0 & p_labels{fi}>0);
        TNs(fi) = sum(a_labels{fi}==0 & p_labels{fi}==0);
        FPs(fi) = sum(a_labels{fi}==0 & p_labels{fi}>0);
        FNs(fi) = sum(a_labels{fi}>0 & p_labels{fi}==0);    
        
        %% QDA --> no hyperparameters to tune
        qda_model = fitcdiscr(x_train, y_train, 'DiscrimType', 'quadratic');

        
        [labels, scores] = predict(qda_model, x_test);
        p_labels_qda{fi} = labels;
        p_scores_qda{fi} = scores;
        a_labels_qda{fi} = y_test;
             
        lengths_qda(fi) = length(p_labels_qda{fi});
        accs_qda(fi) = 100*(1-sum(abs(p_labels_qda{fi}-a_labels_qda{fi}))./length(p_labels_qda{fi}));
        TPs_qda(fi) = sum(a_labels_qda{fi}>0 & p_labels_qda{fi}>0);
        TNs_qda(fi) = sum(a_labels_qda{fi}==0 & p_labels_qda{fi}==0);
        FPs_qda(fi) = sum(a_labels_qda{fi}==0 & p_labels_qda{fi}>0);
        FNs_qda(fi) = sum(a_labels_qda{fi}>0 & p_labels_qda{fi}==0);
        
        toc(tc);
    end
    
    pred_labels = [];
    for fi = 1:no_folds_val
        pred_labels = [pred_labels; p_labels{fi}];
    end
    
    disp(accs)
    disp(sum(accs.*lengths)/sum(lengths))
    disp(std(accs))
    
    disp('best gammas')
    disp(best_gs)
    
    disp('best deltas')
    disp(best_ds)

    
    pred_labels_qda = [];
    for fi = 1:no_folds_val
        pred_labels_qda = [pred_labels_qda; p_labels_qda{fi}];
    end
    
    disp(accs_qda)
    disp(sum(accs_qda.*lengths_qda)/sum(lengths_qda))
    disp(std(accs_qda))
    
