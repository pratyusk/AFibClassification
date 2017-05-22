
function [test_acc, test_sens, test_spec, model] = CART_classifier(X,y,training,test)


   train_x = X(training,:); % form the training data
    train_y = y(training)'; 
    
    test_x = X(test,:); % form the testing data
    test_y = y(test);    
     


    % get RF classifier model
    model = fitctree(train_x,train_y);
    % clculate perfomance measures
         

    
 test_pred_y = model.predict(test_x); % use the trained model to classify the testing data
     

%true pos is class 1 classified as class 1
 true_pos_test = sum(test_y == test_pred_y & test_y == 1);


% true negative is class 0 classified as class 0
true_neg_test = sum(test_y == test_pred_y & test_y == 0);


% false positive is class 0 classified as class 1
false_pos_test = sum(test_y ~= test_pred_y &  test_y == 1);


% false negative is class 1 classified as class 0
false_neg_test = sum(test_y ~= test_pred_y &  test_y == 0);
    
% accuracy = (true pos + true neg)/(true pos + true neg + false pos + false neg)
test_acc = (true_pos_test + true_neg_test)/(true_pos_test + true_neg_test + false_pos_test + false_neg_test);

% sensitivity = true pos/(true pos + false neg) = true pos/ pos class
test_sens = true_pos_test/(true_pos_test + false_neg_test);

% Specificity = true neg/(true neg + false pos)
test_spec = true_neg_test/(true_neg_test + false_pos_test);
    
    
    
    
    
end

