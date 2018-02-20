close all; clear all; clc
%% LOAD DATA
path = ('mixed_cell_pnts_3D_new.mat');
load(path)
fullData = mixed_cell_pnts_new;

%% CALL METHOD
%% INPUT: DATA (2D xy coordinate)
%         plotflag: 1 is for ploting, 0 is for non-ploting
%  OUTPUT: centre (2D xy value)

% M_final = FindCenter(fullData,1);
for i = 1: length(fullData)
    cell = fullData{i}(:,1:2);
%     figure; scatter(cell(:,1),cell(:,2));
    flag = checkCellisCircle(cell);
    if(flag == true)
        fprintf('cell %i is circle cell \n',i)        
    else
        fprintf('cell %i is not circle cell \n',i)        
    end
end