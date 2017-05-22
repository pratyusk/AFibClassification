%% nested cross validation to train svm classifier 

    clear all; close all; clc;
    no_folds_val = 10;
    load('Phi_train_new.mat');  
    features = [2 3 5 7 8 10];
    data = Phi(:,[features+1, 1]);
    data(find(data(:,end) == 65),end) = 1;
    data(find(data(:,end) == 78),end) = 0;    

    % Hyperparameter
    NumTrees = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 25 30 35 40 45 50];
    
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
    best_num_trees = zeros(1,no_folds_val);
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
        
        acc_fold = zeros(length(NumTrees),1);
        
        % inner cross-validation
        inner_data0 = [x_train(y_train==0,:) y_train(y_train==0,:)];
        inner_data1 = [x_train(y_train==1,:) y_train(y_train==1,:)];
    
        
        inner_rperm0 = randperm(length(inner_data0));
        inner_rperm1 = randperm(length(inner_data1));
        
        inner_data0 = inner_data0(inner_rperm0,:);
        inner_data1 = inner_data1(inner_rperm1,:);

        inner_fold_size0 = ceil(size(inner_data0,1)/no_folds_val);
        inner_fold_size1 = ceil(size(inner_data1,1)/no_folds_val); 
        inner_accs = zeros(1,10);
        for ti = 1:length(NumTrees)
            for k = 1:no_folds_val
                inner_fold_test0 = inner_data0(inner_fold_size0*(k-1)+1:min(inner_fold_size0*k,size(inner_data0,1)),:);
                inner_fold_test1 = inner_data1(inner_fold_size1*(k-1)+1:min(inner_fold_size1*k,size(inner_data1,1)),:);

                inner_fold_train0 = inner_data0([1:inner_fold_size0*(k-1),min(inner_fold_size0*k,size(inner_data0,1))+1:end],:);
                inner_fold_train1 = inner_data1([1:inner_fold_size1*(k-1),min(inner_fold_size1*k,size(inner_data1,1))+1:end],:);            

                inner_fold_test = [inner_fold_test0 ; inner_fold_test1];
                inner_fold_test = inner_fold_test(randperm(size(inner_fold_test,1)),:);

                inner_fold_train = [inner_fold_train0 ; inner_fold_train1];
                inner_fold_train = inner_fold_train(randperm(size(inner_fold_train,1)),:);

                inner_x_train = inner_fold_train(:,1:end-1);
                inner_y_train = inner_fold_train(:,end);
                inner_x_test = inner_fold_test(:,1:end-1);
                inner_y_test = inner_fold_test(:,end);
                tree_model = TreeBagger(NumTrees(ti),inner_x_train, inner_y_train);
                [inner_labels, ~] = tree_model.predict(inner_x_test);
                inner_accs(k) = 100*(1-sum(abs(str2num(cell2mat(inner_labels))-inner_y_test))./length(inner_labels));
            end
            acc_fold(ti) = mean(inner_accs);
        end
        [mt,mi] = min(acc_fold);
        best_num_tree = NumTrees(mi);
        
        tree_model = TreeBagger(best_num_tree,x_train, y_train);
        
        [labels, scores] = tree_model.predict(x_test);
        p_labels{fi} = str2num(cell2mat(labels));
        p_scores{fi} = scores;
        a_labels{fi} = y_test;
        best_num_trees(fi) = best_num_tree;
       
        
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
    
%     save('resutls','p_labels','p_scores','a_labels','best_cs','best_ss');
