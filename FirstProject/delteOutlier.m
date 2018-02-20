function XX = delteOutlier(X)
    X = unique(X,'rows');
    D = pdist(X);
    Z = squareform(D);
%     maxvector = (max(Z));
%     meanvector = (mean(Z));
    for i=1:size(Z,2)
        Z(i,i)=9999;
    end
    minvector = (min(Z));
    idx = isoutlier(minvector);
    XX = X; XX(idx==1,:) = [];
    
%     idx = isoutlier(maxvector);
%     XX = X; XX(idx==1,:) = [];
%     
%     idx = isoutlier(meanvector);
%     XX = X; XX(idx==1,:) = [];

end