function svm_plot_classes(X,y,classifier)

    pred = predict(classifier,X); % find the predictions
        
    [xGrid,yGrid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),100),linspace(min(X(:,2)),max(X(:,2)),100)); 
    % create a mesh grid with x and y values covering the surface of the
    % plot
    
    [~,scores] = predict(classifier,[xGrid(:),yGrid(:)]); % find the classifier output for all the points in the mesh grid  
    scores = reshape(scores(:,1),size(xGrid)); % reshape the classifier output back into a matrix form
    
    t = (pred == y); % find the correctly classified data points
    f = (pred ~= y); % find the incorrectly classified data points
    
    figure;        
    contourf(xGrid,yGrid,scores,[0 0]); % draw the contour showing the class boundary
    
    hold all;
    h1 = gscatter(X(t,1),X(t,2),y(t)); % plot the correctly classified data points
    hold all;
    h2 = gscatter(X(f,1),X(f,2),y(f),[],'**'); % plot the incorrectly classified data points
    lgnds = {'Class 1 - Correctly Classified', 'Class 2 - Correctly Classified', 'Class 1 - Incorrectly Classified', 'Class 2 - Incorrectly Classified'};
    lgnds = lgnds([any(y(t)==1), any(y(t)==2), any(y(f)==1), any(y(f)==2)]); % select legend entries for which there is at least one point in the plot
    legend([h1;h2],lgnds);
    axis tight
    
    colormap(autumn); % change the color of the contour
    caxis([-10,0.01])

end

