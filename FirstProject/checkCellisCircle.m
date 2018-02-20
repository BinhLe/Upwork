function flag = checkCellisCircle(cell)
cell = delteOutlier(unique(cell,'rows'));

if (length(cell)<50)
    flag = false;
else
    [f,xi,bw] = ksdensity(cell);
    meanbw = mean(bw);
    % disp(meanbw);
    if (meanbw >=0.0010)
        flag = true;
    else
        flag = false;
    end
end