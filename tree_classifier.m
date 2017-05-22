function [test_acc, train_acc, model] = tree_classifier(X,y,max_n_splits,training,test)
    
    train_X = X(training,:); % form the training data
    train_y = y(training); 
    
    test_X = X(test,:); % form the testing data
    test_y = y(test);         

    model = fitctree(train_X, train_y, 'MaxNumSplits', max_n_splits);
    % train a decision tree model with a non-linear kernel and its kernel parameter using the training data

    train_pred_y = predict(model, train_X); % use the trained model to classify the training data
    train_acc = sum(train_y == train_pred_y)/length(train_y); % find the training accuracy               
    
    test_pred_y = predict(model, test_X); % use the trained model to classify the testing data
    test_acc = sum(test_y == test_pred_y)/length(test_y); % find the test accuracy
     
end

