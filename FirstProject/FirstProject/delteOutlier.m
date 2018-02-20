function XX = delteOutlier(X)
    D = pdist(X);
    Z = squareform(D);
    for i=1:size(Z,2)
        Z(i,i)=9999;
    end
    minvector = (min(Z));
    idx = isoutlier(minvector);
    XX = X; XX(idx==1,:) = [];

end